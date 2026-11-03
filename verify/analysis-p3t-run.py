"""Run a command, report wall time, peak working set (Windows) and exit code.  python -I verify/analysis-p3t-run.py LABEL STDOUT_FILE -- cmd args...
Prints one line:  RUN | LABEL | rc=<n> | wall=<s> s | peak_ws=<MB> MB"""
import ctypes, ctypes.wintypes as W, subprocess, sys, time

class PMC(ctypes.Structure):
    _fields_ = [("cb", W.DWORD), ("PageFaultCount", W.DWORD), ("PeakWorkingSetSize", ctypes.c_size_t), ("WorkingSetSize", ctypes.c_size_t),
                ("QuotaPeakPagedPoolUsage", ctypes.c_size_t), ("QuotaPagedPoolUsage", ctypes.c_size_t), ("QuotaPeakNonPagedPoolUsage", ctypes.c_size_t),
                ("QuotaNonPagedPoolUsage", ctypes.c_size_t), ("PagefileUsage", ctypes.c_size_t), ("PeakPagefileUsage", ctypes.c_size_t)]

label, outf = sys.argv[1], sys.argv[2]
cmd = sys.argv[sys.argv.index("--") + 1:]
t0 = time.time()
with open(outf, "wb") as fo:
    p = subprocess.Popen(cmd, stdout=fo, stderr=subprocess.STDOUT)
    rc = p.wait()
wall = time.time() - t0
pm = PMC(); pm.cb = ctypes.sizeof(pm)
k = ctypes.WinDLL("kernel32", use_last_error=True); ps = ctypes.WinDLL("psapi", use_last_error=True)
k.GetCurrentProcess.restype = W.HANDLE
ps.GetProcessMemoryInfo.argtypes = [W.HANDLE, ctypes.POINTER(PMC), W.DWORD]
ok = ps.GetProcessMemoryInfo(W.HANDLE(int(p._handle)), ctypes.byref(pm), pm.cb)
peak = pm.PeakWorkingSetSize / 2**20 if ok else float("nan")
print("RUN | %s | rc=%d | wall=%.1f s | peak_ws=%.0f MB" % (label, rc, wall, peak))
