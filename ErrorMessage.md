[100%] Linking CXX executable zhijian
/usr/bin/ld: /home/z651w035/hpmusic/contrib/cgnslib/lib/libcgns.a(ADFH.c.o): in function `ADFH_Children_IDs':
ADFH.c:(.text+0x441f): undefined reference to `H5Literate'
/usr/bin/ld: ADFH.c:(.text+0x453b): undefined reference to `H5Literate'
/usr/bin/ld: /home/z651w035/hpmusic/contrib/cgnslib/lib/libcgns.a(ADFH.c.o): in function `ADFH_Children_Names':
ADFH.c:(.text+0x581f): undefined reference to `H5Literate'
/usr/bin/ld: ADFH.c:(.text+0x592b): undefined reference to `H5Literate'
collect2: error: ld returned 1 exit status
make[2]: *** [CMakeFiles/zhijian.dir/build.make:157: zhijian] Error 1
make[1]: *** [CMakeFiles/Makefile2:111: CMakeFiles/zhijian.dir/all] Error 2
make: *** [Makefile:136: all] Error 2

