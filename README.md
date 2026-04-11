# nanoRQc

nanoRQc is a simpler version of nanorq for the convenience. 
It only supports core feature (encoder and decoder) of RaptorQ.

While nanorq is compact and performant, the API isn't easy to use. 
So, I refined the library for my usage. 
Thanks Joseph Calderon for releasing the original nanorq library.

# Usable on Microsoft Visual Studio

This is a modified version of [nanorq](https://github.com/sleepybishop/nanorq). 
Because the original nanorq cannot be compiled on Microsoft Visual Studio, 
I changed some points for compatibility with Microsoft Visual C.

I removed mmap feature, as C-runtime library doesn't support it on Windows OS.

I added detecting CPU feature (SSSE3 or AVX2) to switch functions automatically.

# Convenient API for utility

I removed block selection. Source Block Number is 0 always. 
To support large data in single block, symbol size can be larger than 65536 bytes. 
Only two parameters (K and T) are used to setup encoder.

I removed support of file stream. 
Because it only supports data on memory, decoding may become faster.

I added more return values to determine failure. 
It may be helpful to control the behavior after error.

I added some useful functions. 
It's possible to add each source symbol in encoder one by one, or add all source symbols at once. 
It's possible to create each repair symbol one by one, or create all repair symbols at once. 
It's possible to restore each symbol in decoder one by one, or recover all missing source symbols at once. 
You may try the possibility of successful repair, while adding symbols into decoder.

# Notes

As I changed some, gcc may not compile this modified version. 
If you use nanorq with gcc on Linux or Unix, you would better use [original version](https://github.com/sleepybishop/nanorq).

If you want full support of RFC 6330 on Microsoft Visual Studio, 
you may use [nanorq for Visual C](https://github.com/Yutaka-Sawada/nanorq_vc).

