#!/usr/bin/env python3
"""Build and install REBOUNDx 5.1.0 and celmech 1.5.8 on native Windows (x64, CPython 3.12).

Starting from the pristine PyPI sdists this script
  1. creates (or reuses) a uv virtual environment and installs the build requirements,
  2. downloads the sdists from PyPI (sha256 pinned) and extracts them into a fresh work directory,
  3. applies patches/reboundx-5.1.0-windows.diff, patches/celmech-1.5.8-windows.diff and (unless --no-compat)
     patches/celmech-1.5.8-modern-deps-compat.diff,
  4. builds both wheels (REBOUNDx with MSVC, celmech's libcelmech with `python -m ziglang cc`),
  5. installs the wheels (plus celmech's dependencies) into the venv.

Requirements: Windows x64, `uv` on PATH, `git` on PATH (only for `git apply`), MSVC Build Tools 2022
(found automatically by setuptools). Nothing is installed outside the venv and the work directory.
Rerunnable: the source trees are re-extracted and re-patched on every run.

    python build_windows.py [--venv DIR] [--work DIR] [--wheels DIR] [--python 3.12]
"""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tarfile
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))

# name -> (version, sha256 of the PyPI sdist)
SDISTS = {
    'reboundx': ('5.1.0', '25fcccc8194a83b548a8e11d81f295fb12a21ece2538e61be93d8256d9cbc1b2'),
    'celmech': ('1.5.8', 'b40743195571fc39d187b2d41eca5357c8be181dd8d74c2e1dce4360fe48f399'),
}
# Build requirements (setup.py of both packages imports rebound/numpy; celmech's libcelmech needs ziglang)
# plus `pytest` (to run celmech's tests) and `ipython`, which celmech imports (celmech/poincare.py) but does not declare as a dependency.
BUILD_REQS = ['rebound>=5.0.0', 'numpy', 'scipy', 'sympy', 'mpmath', 'setuptools>=68', 'wheel', 'ziglang',
              'ipython', 'pytest']


def run(cmd, **kw):
    print('+', ' '.join(map(str, cmd)), flush=True)
    subprocess.check_call(cmd, **kw)


def fetch_sdist(name, version, sha256, dl_dir):
    path = os.path.join(dl_dir, f'{name}-{version}.tar.gz')
    if not os.path.exists(path):
        with urllib.request.urlopen(f'https://pypi.org/pypi/{name}/{version}/json') as r:
            info = json.load(r)
        url = next(u['url'] for u in info['urls'] if u['packagetype'] == 'sdist')
        print('downloading', url, flush=True)
        with urllib.request.urlopen(url) as r, open(path, 'wb') as f:
            shutil.copyfileobj(r, f)
    digest = hashlib.sha256(open(path, 'rb').read()).hexdigest()
    if digest != sha256:
        raise SystemExit(f'sha256 mismatch for {path}: {digest} != {sha256}')
    return path


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--venv', default=os.path.join(HERE, 'venv'))
    ap.add_argument('--work', default=os.path.join(HERE, '_work'))
    ap.add_argument('--wheels', default=os.path.join(HERE, 'wheels'))
    ap.add_argument('--python', default='3.12', help='interpreter request used when the venv must be created')
    ap.add_argument('--no-compat', action='store_true',
                    help='skip celmech-1.5.8-modern-deps-compat.diff (REBOUND 5 / SciPy / Python 3 fixes unrelated '
                         'to Windows)')
    args = ap.parse_args()
    if sys.platform != 'win32':
        raise SystemExit('this script is for native Windows only')
    uv = shutil.which('uv') or sys.exit('uv not found on PATH')
    git = shutil.which('git') or sys.exit('git not found on PATH (needed for `git apply`)')
    venv, work, wheels = (os.path.abspath(p) for p in (args.venv, args.work, args.wheels))
    os.makedirs(work, exist_ok=True)
    os.makedirs(wheels, exist_ok=True)

    if not os.path.exists(os.path.join(venv, 'Scripts', 'python.exe')):
        run([uv, 'venv', '-p', args.python, venv])
    py = os.path.join(venv, 'Scripts', 'python.exe')
    run([uv, 'pip', 'install', '-p', py] + BUILD_REQS)

    # Environment for builds and patching:
    #  - GIT_CEILING_DIRECTORIES keeps git from finding an unrelated repository above the work dir (otherwise
    #    setup.py would embed that repository's HEAD as the library git hash, and `git apply` would misbehave).
    #  - zig keeps its caches inside the work dir instead of %LOCALAPPDATA%.
    env = dict(os.environ, GIT_CEILING_DIRECTORIES=work,
               ZIG_GLOBAL_CACHE_DIR=os.path.join(work, 'zigcache', 'global'),
               ZIG_LOCAL_CACHE_DIR=os.path.join(work, 'zigcache', 'local'))

    dl_dir = os.path.join(work, 'dl')
    src_root = os.path.join(work, 'src')
    os.makedirs(dl_dir, exist_ok=True)
    if os.path.isdir(src_root):
        shutil.rmtree(src_root)
    os.makedirs(src_root)

    # reboundx first: celmech requires reboundx>=4.0.0 and must not pull the (unbuildable) PyPI sdist.
    for name in ('reboundx', 'celmech'):
        version, sha = SDISTS[name]
        tarball = fetch_sdist(name, version, sha, dl_dir)
        with tarfile.open(tarball) as t:
            t.extractall(src_root, filter='data')
        src = os.path.join(src_root, f'{name}-{version}')
        patches = [f'{name}-{version}-windows.diff']
        if name == 'celmech' and not args.no_compat:
            patches.append(f'{name}-{version}-modern-deps-compat.diff')
        for pf in patches:
            run([git, 'apply', '--whitespace=nowarn', '-p1', os.path.join(HERE, 'patches', pf)], cwd=src, env=env)
        run([uv, 'build', '--wheel', '--no-build-isolation', '-p', py, '--out-dir', wheels, src], env=env)

    whl = {n: os.path.join(wheels, f'{n}-{v}-cp312-cp312-win_amd64.whl') for n, (v, _) in SDISTS.items()}
    for w in whl.values():
        if not os.path.exists(w):
            raise SystemExit(f'expected wheel missing: {w}')
    run([uv, 'pip', 'install', '-p', py, '--reinstall-package', 'reboundx', '--reinstall-package', 'celmech',
         whl['reboundx'], whl['celmech'], 'ipython'], env=env)
    print('\nDone. Test with:\n  ' + py + ' -c "import rebound, reboundx, celmech; '
          'print(rebound.__version__, reboundx.__version__, celmech.__version__)"')


if __name__ == '__main__':
    main()
