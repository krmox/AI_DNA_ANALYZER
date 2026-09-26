"""V2 BLAKE2b key == hashlib.blake2b(digest_size=8) big-endian, incl. block-boundary lengths."""
import ctypes, hashlib, os, subprocess, sys, tempfile, textwrap
src = textwrap.dedent('''
#include <cstdio>
#include <cstring>
#include "dnav2/blake2b.hpp"
extern "C" unsigned long long key(const char* s, unsigned long n){ return dnav2::blake2b_key64((const uint8_t*)s, n); }
''')
d = tempfile.mkdtemp(dir="/mnt/archive/AI_DNA_ANALYZER_v2/work")
open(f"{d}/b.cpp", "w").write(src)
inc = os.path.join(os.path.dirname(__file__), "..", "include")
subprocess.check_call(["g++", "-O2", "-std=c++20", "-shared", "-fPIC", f"-I{inc}", "-o", f"{d}/b.so", f"{d}/b.cpp"])
lib = ctypes.CDLL(f"{d}/b.so"); lib.key.restype = ctypes.c_ulonglong; lib.key.argtypes = [ctypes.c_char_p, ctypes.c_ulong]
import random
random.seed(1)
bad = 0; n = 0
names = [b"", b"a", b"x" * 127, b"x" * 128, b"x" * 129, b"x" * 255, b"x" * 256, b"x" * 257]
names += [bytes(random.choices(b"ABCDEFGHIJ0123456789:_-", k=random.randint(1, 300))) for _ in range(20000)]
for m in names:
    ref = int.from_bytes(hashlib.blake2b(m, digest_size=8).digest(), "big")
    n += 1; bad += lib.key(m, len(m)) != ref
print(f"[{'PASS' if bad == 0 else 'FAIL'}] blake2b_key_equals_hashlib: n={n} mismatches={bad}")
if os.environ.get("V2_CORRECTNESS_CSV"):
    import time
    open(os.environ["V2_CORRECTNESS_CSV"], "a").write(f"kernel_tests,{time.strftime('%F')},blake2b_key_equals_hashlib,{n},{bad},{'PASS' if bad==0 else 'FAIL'},\n")
sys.exit(1 if bad else 0)
