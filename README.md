# nanoRQc

nanoRQc is a simpler version of nanorq for the convenience. 
It only supports core feature (encoder and decoder) of RaptorQ.

While nanorq is compact and performant, the API isn't easy to use. 
So, I refined the library for my usage. 
Thanks Joseph Calderon for releasing the original nanorq library.

# Usable on Microsoft Visual Studio

This is a modified version of [nanorq](https://github.com/ursoid-tech/nanorq). 
Because the original nanorq cannot be compiled on Microsoft Visual Studio, 
I changed some points for compatibility with Microsoft Visual C.

To enable atomic functions, you must use Microsoft Visual Studio 2022 (version 17.5 or later).  
It requires compiler option "`/std:c11`" or "`/std:c17`".  
It requires "`/experimental:c11atomics`" as additional compiler option, too.

# Convenient API for utility

I removed block selection. Source Block Number is 0 always. 
To support large data in single block, symbol size can be larger than 65536 bytes. 
Only two parameters (K and T) are used to setup encoder.

I removed support of file stream and memory mapped file. 

I added more return values to determine failure. 
It may be helpful to control the behavior after error.

I added some useful functions. 
It's possible to add each source symbol in encoder one by one. 
It's possible to restore each symbol in decoder one by one. 
You may try the possibility of successful repair, while adding symbols into decoder.

# Notes

As I changed some, gcc may not compile this modified version. 
If you use nanorq with gcc on Linux or Unix, you would better use [original version](https://github.com/ursoid-tech/nanorq).

If you want full support of RFC 6330 on Microsoft Visual Studio, 
you may use [nanorq for Visual C](https://github.com/Yutaka-Sawada/nanorq_vc).

