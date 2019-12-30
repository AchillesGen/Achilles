import ctypes

libc = ctypes.CDLL('libc.so.6')

print(libc.time(None))
