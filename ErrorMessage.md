[ 92%] Linking CUDA device code CMakeFiles/zhijian.dir/cmake_device_link.o
nvlink warning : Skipping incompatible '/usr/lib64/libdl.a' when searching for -ldl (target: sm_70)
nvlink warning : Skipping incompatible '/usr/lib64/librt.a' when searching for -lrt (target: sm_70)
nvlink warning : Skipping incompatible '/usr/lib64/libpthread.a' when searching for -lpthread (target: sm_70)
nvlink warning : Skipping incompatible '/usr/lib64/libdl.a' when searching for -ldl (target: sm_75)
nvlink warning : Skipping incompatible '/usr/lib64/librt.a' when searching for -lrt (target: sm_75)
nvlink warning : Skipping incompatible '/usr/lib64/libpthread.a' when searching for -lpthread (target: sm_75)
nvlink warning : Skipping incompatible '/usr/lib64/libdl.a' when searching for -ldl (target: sm_80)
nvlink warning : Skipping incompatible '/usr/lib64/librt.a' when searching for -lrt (target: sm_80)
nvlink warning : Skipping incompatible '/usr/lib64/libpthread.a' when searching for -lpthread (target: sm_80)
nvlink warning : Skipping incompatible '/usr/lib64/libdl.a' when searching for -ldl (target: sm_86)
nvlink warning : Skipping incompatible '/usr/lib64/librt.a' when searching for -lrt (target: sm_86)
nvlink warning : Skipping incompatible '/usr/lib64/libpthread.a' when searching for -lpthread (target: sm_86)
nvlink warning : Skipping incompatible '/usr/lib64/libdl.a' when searching for -ldl (target: sm_89)
nvlink warning : Skipping incompatible '/usr/lib64/librt.a' when searching for -lrt (target: sm_89)
nvlink warning : Skipping incompatible '/usr/lib64/libpthread.a' when searching for -lpthread (target: sm_89)
nvlink warning : Skipping incompatible '/usr/lib64/libdl.a' when searching for -ldl (target: sm_90)
nvlink warning : Skipping incompatible '/usr/lib64/librt.a' when searching for -lrt (target: sm_90)
nvlink warning : Skipping incompatible '/usr/lib64/libpthread.a' when searching for -lpthread (target: sm_90)
[100%] Linking CXX executable zhijian
/usr/bin/ld: /home/z651w035/hpmusic/contrib/cgnslib/lib/libcgns.a(ADFH.c.o): in function `ADFH_Children_IDs':
ADFH.c:(.text+0x441f): undefined reference to `H5Literate'
/usr/bin/ld: ADFH.c:(.text+0x453b): undefined reference to `H5Literate'
/usr/bin/ld: /home/z651w035/hpmusic/contrib/cgnslib/lib/libcgns.a(ADFH.c.o): in function `ADFH_Children_Names':
ADFH.c:(.text+0x581f): undefined reference to `H5Literate'
/usr/bin/ld: ADFH.c:(.text+0x592b): undefined reference to `H5Literate'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__ismalloc':
gklib.c:(.text+0x10c5): undefined reference to `gk_malloc'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__iAllocMatrix':
gklib.c:(.text+0x1142): undefined reference to `gk_malloc'
/usr/bin/ld: gklib.c:(.text+0x119c): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__iFreeMatrix':
gklib.c:(.text+0x11eb): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__rsmalloc':
gklib.c:(.text+0x1333): undefined reference to `gk_malloc'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__rAllocMatrix':
gklib.c:(.text+0x13a4): undefined reference to `gk_malloc'
/usr/bin/ld: gklib.c:(.text+0x13fc): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__rFreeMatrix':
gklib.c:(.text+0x144b): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__ikvsmalloc':
gklib.c:(.text+0x1515): undefined reference to `gk_malloc'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__ikvAllocMatrix':
gklib.c:(.text+0x1592): undefined reference to `gk_malloc'
/usr/bin/ld: gklib.c:(.text+0x15ec): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__ikvFreeMatrix':
gklib.c:(.text+0x163b): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__rkvsmalloc':
gklib.c:(.text+0x1705): undefined reference to `gk_malloc'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__rkvAllocMatrix':
gklib.c:(.text+0x1782): undefined reference to `gk_malloc'
/usr/bin/ld: gklib.c:(.text+0x17dc): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__rkvFreeMatrix':
gklib.c:(.text+0x182b): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__ipqInit':
gklib.c:(.text+0x18dc): undefined reference to `gk_idxsmalloc'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__ipqCreate':
gklib.c:(.text+0x1907): undefined reference to `gk_malloc'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__ipqFree':
gklib.c:(.text+0x1976): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__ipqDestroy':
gklib.c:(.text+0x19ad): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__rpqInit':
gklib.c:(.text+0x1e9c): undefined reference to `gk_idxsmalloc'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__rpqCreate':
gklib.c:(.text+0x1ec7): undefined reference to `gk_malloc'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__rpqFree':
gklib.c:(.text+0x1f36): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__rpqDestroy':
gklib.c:(.text+0x1f6d): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__iargmax_n':
gklib.c:(.text+0x38ae): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__rargmax_n':
gklib.c:(.text+0x3c52): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__imalloc':
gklib.c:(.text+0xfc5): undefined reference to `gk_malloc'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__irealloc':
gklib.c:(.text+0xfd5): undefined reference to `gk_realloc'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__iFreeMatrix':
gklib.c:(.text+0x120a): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__rmalloc':
gklib.c:(.text+0x1255): undefined reference to `gk_malloc'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__rrealloc':
gklib.c:(.text+0x1265): undefined reference to `gk_realloc'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__rFreeMatrix':
gklib.c:(.text+0x146a): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__ikvmalloc':
gklib.c:(.text+0x14c5): undefined reference to `gk_malloc'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__ikvrealloc':
gklib.c:(.text+0x14d5): undefined reference to `gk_realloc'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__ikvFreeMatrix':
gklib.c:(.text+0x165a): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__rkvmalloc':
gklib.c:(.text+0x16b5): undefined reference to `gk_malloc'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__rkvrealloc':
gklib.c:(.text+0x16c5): undefined reference to `gk_realloc'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__rkvFreeMatrix':
gklib.c:(.text+0x184a): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__isrand':
gklib.c:(.text+0x2454): undefined reference to `gk_randinit'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(gklib.c.o): in function `libmetis__irand':
gklib.c:(.text+0x2461): undefined reference to `gk_randint32'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(kmetis.c.o): in function `libmetis__InitKWayPartitioning':
kmetis.c:(.text+0xf0): undefined reference to `gk_errexit'
/usr/bin/ld: kmetis.c:(.text+0x103): undefined reference to `gk_free'
/usr/bin/ld: kmetis.c:(.text+0x182): undefined reference to `gk_errexit'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(kmetis.c.o): in function `libmetis__MlevelKWayPartitioning':
kmetis.c:(.text+0x1df): undefined reference to `gk_errexit'
/usr/bin/ld: kmetis.c:(.text+0x339): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: kmetis.c:(.text+0x359): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(kmetis.c.o): in function `METIS_PartGraphKway':
kmetis.c:(.text+0x1401): undefined reference to `gk_malloc_init'
/usr/bin/ld: kmetis.c:(.text+0x1410): undefined reference to `gk_sigtrap'
/usr/bin/ld: kmetis.c:(.text+0x1418): undefined reference to `gk_cur_jbufs'
/usr/bin/ld: kmetis.c:(.text+0x142b): undefined reference to `gk_jbufs'
/usr/bin/ld: kmetis.c:(.text+0x1458): undefined reference to `gk_siguntrap'
/usr/bin/ld: kmetis.c:(.text+0x145f): undefined reference to `gk_malloc_cleanup'
/usr/bin/ld: kmetis.c:(.text+0x1554): undefined reference to `gk_log2'
/usr/bin/ld: kmetis.c:(.text+0x1668): undefined reference to `gk_log2'
/usr/bin/ld: kmetis.c:(.text+0x1695): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: kmetis.c:(.text+0x16c1): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: kmetis.c:(.text+0x16f3): undefined reference to `gk_siguntrap'
/usr/bin/ld: kmetis.c:(.text+0x1717): undefined reference to `gk_errexit'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(kwayrefine.c.o): in function `libmetis__AllocateKWayPartitionMemory':
kwayrefine.c:(.text+0x9f): undefined reference to `gk_malloc'
/usr/bin/ld: kwayrefine.c:(.text+0xcb): undefined reference to `gk_malloc'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(kwayrefine.c.o): in function `libmetis__ComputeKWayBoundary':
kwayrefine.c:(.text+0x144): undefined reference to `gk_errexit'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(kwayrefine.c.o): in function `libmetis__ProjectKWayPartition':
kwayrefine.c:(.text+0x134f): undefined reference to `gk_errexit'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(kwayrefine.c.o): in function `libmetis__RefineKWay':
kwayrefine.c:(.text+0x17d1): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: kwayrefine.c:(.text+0x17f1): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: kwayrefine.c:(.text+0x1871): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: kwayrefine.c:(.text+0x1891): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: kwayrefine.c:(.text+0x19b1): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(kwayrefine.c.o):kwayrefine.c:(.text+0x1b11): more undefined references to `gk_CPUSeconds' follow
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(kwayrefine.c.o): in function `libmetis__AllocateKWayPartitionMemory':
kwayrefine.c:(.text+0x89): undefined reference to `gk_errexit'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(kwayrefine.c.o): in function `libmetis__ComputeKWayPartitionParams':
kwayrefine.c:(.text+0xa0f): undefined reference to `gk_errexit'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(minconn.c.o): in function `libmetis__ComputeSubDomainGraph':
minconn.c:(.text+0xe9): undefined reference to `gk_errexit'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(minconn.c.o): in function `libmetis__EliminateSubDomainEdges':
minconn.c:(.text+0x1b39): undefined reference to `gk_errexit'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(minconn.c.o): in function `libmetis__PrintSubDomainGraph':
minconn.c:(.text+0x2380): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(options.c.o): in function `libmetis__FreeCtrl':
options.c:(.text+0xb38): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(options.c.o): in function `libmetis__SetupCtrl':
options.c:(.text+0xb8d): undefined reference to `gk_malloc'
/usr/bin/ld: options.c:(.text+0xbe5): undefined reference to `gk_errexit'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(pmetis.c.o): in function `libmetis__SplitGraphPart':
pmetis.c:(.text+0xae9): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: pmetis.c:(.text+0xb08): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(pmetis.c.o): in function `METIS_PartGraphRecursive':
pmetis.c:(.text+0xeb1): undefined reference to `gk_malloc_init'
/usr/bin/ld: pmetis.c:(.text+0xec0): undefined reference to `gk_sigtrap'
/usr/bin/ld: pmetis.c:(.text+0xec8): undefined reference to `gk_cur_jbufs'
/usr/bin/ld: pmetis.c:(.text+0xedb): undefined reference to `gk_jbufs'
/usr/bin/ld: pmetis.c:(.text+0xf08): undefined reference to `gk_siguntrap'
/usr/bin/ld: pmetis.c:(.text+0xf0f): undefined reference to `gk_malloc_cleanup'
/usr/bin/ld: pmetis.c:(.text+0x1081): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: pmetis.c:(.text+0x10a9): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: pmetis.c:(.text+0x10d3): undefined reference to `gk_siguntrap'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(refine.c.o): in function `libmetis__Refine2Way':
refine.c:(.text+0xc43): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: refine.c:(.text+0xcb0): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: refine.c:(.text+0xcd1): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: refine.c:(.text+0xd01): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: refine.c:(.text+0xd21): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(refine.c.o):refine.c:(.text+0xd41): more undefined references to `gk_CPUSeconds' follow
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(wspace.c.o): in function `libmetis__AllocateWorkSpace':
wspace.c:(.text+0x2f): undefined reference to `gk_mcoreCreate'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(wspace.c.o): in function `libmetis__AllocateRefinementWorkSpace':
wspace.c:(.text+0xca): undefined reference to `gk_errexit'
/usr/bin/ld: wspace.c:(.text+0x19c): undefined reference to `gk_malloc'
/usr/bin/ld: wspace.c:(.text+0x1c0): undefined reference to `gk_malloc'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(wspace.c.o): in function `libmetis__FreeWorkSpace':
wspace.c:(.text+0x1e5): undefined reference to `gk_mcoreDestroy'
/usr/bin/ld: wspace.c:(.text+0x202): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(wspace.c.o): in function `libmetis__cnbrpoolGetNext':
wspace.c:(.text+0x39c): undefined reference to `gk_realloc'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(wspace.c.o): in function `libmetis__vnbrpoolGetNext':
wspace.c:(.text+0x44f): undefined reference to `gk_realloc'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(wspace.c.o): in function `libmetis__FreeWorkSpace':
wspace.c:(.text+0x2a5): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(wspace.c.o): in function `libmetis__wspacemalloc':
wspace.c:(.text+0x2b8): undefined reference to `gk_mcoreMalloc'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(wspace.c.o): in function `libmetis__wspacepush':
wspace.c:(.text+0x2c8): undefined reference to `gk_mcorePush'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(wspace.c.o): in function `libmetis__wspacepop':
wspace.c:(.text+0x2d8): undefined reference to `gk_mcorePop'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(coarsen.c.o): in function `libmetis__Match_2HopAny':
coarsen.c:(.text+0x341): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: coarsen.c:(.text+0x361): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(coarsen.c.o): in function `libmetis__Match_2HopAll':
coarsen.c:(.text+0x708): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: coarsen.c:(.text+0x724): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(coarsen.c.o): in function `libmetis__CreateCoarseGraph':
coarsen.c:(.text+0x18d2): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(coarsen.c.o):coarsen.c:(.text+0x18f5): more undefined references to `gk_CPUSeconds' follow
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(coarsen.c.o): in function `libmetis__CoarsenGraph':
coarsen.c:(.text+0x28fb): undefined reference to `gk_errexit'
/usr/bin/ld: coarsen.c:(.text+0x2a11): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: coarsen.c:(.text+0x2a31): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(coarsen.c.o): in function `CoarsenGraphNlevels':
coarsen.c:(.text+0x2b3b): undefined reference to `gk_errexit'
/usr/bin/ld: coarsen.c:(.text+0x2c61): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: coarsen.c:(.text+0x2c81): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(coarsen.c.o): in function `libmetis__Match_JC':
coarsen.c:(.text+0x3189): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: coarsen.c:(.text+0x31b6): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(contig.c.o): in function `libmetis__FindPartitionInducedComponents':
contig.c:(.text+0x1f4): undefined reference to `gk_free'
/usr/bin/ld: contig.c:(.text+0x219): undefined reference to `gk_free'
/usr/bin/ld: contig.c:(.text+0x22f): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(contig.c.o): in function `libmetis__IsConnectedSubdomain':
contig.c:(.text+0x8ed): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(contig.c.o): in function `libmetis__FindSepInducedComponents':
contig.c:(.text+0xd04): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(contig.c.o): in function `libmetis__EliminateComponents':
contig.c:(.text+0x1f4c): undefined reference to `gk_errexit'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(debug.c.o): in function `libmetis__ComputeVolume':
debug.c:(.text+0x4ca): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(debug.c.o): in function `libmetis__ComputeMaxCut':
debug.c:(.text+0x5e1): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(graph.c.o): in function `libmetis__CreateGraph':
graph.c:(.text+0x1ef): undefined reference to `gk_malloc'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(graph.c.o): in function `libmetis__FreeSData':
graph.c:(.text+0x619): undefined reference to `gk_free'
/usr/bin/ld: graph.c:(.text+0x62d): undefined reference to `gk_free'
/usr/bin/ld: graph.c:(.text+0x641): undefined reference to `gk_free'
/usr/bin/ld: graph.c:(.text+0x655): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(graph.c.o): in function `libmetis__FreeRData':
graph.c:(.text+0x6cb): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(graph.c.o):graph.c:(.text+0x72f): more undefined references to `gk_free' follow
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(graph.c.o): in function `libmetis__graph_WriteToDisk':
graph.c:(.text+0x7b9): undefined reference to `gk_rmpath'
/usr/bin/ld: graph.c:(.text+0x905): undefined reference to `gk_rmpath'
/usr/bin/ld: graph.c:(.text+0x9c9): undefined reference to `gk_free'
/usr/bin/ld: graph.c:(.text+0x9e1): undefined reference to `gk_free'
/usr/bin/ld: graph.c:(.text+0x9f9): undefined reference to `gk_free'
/usr/bin/ld: graph.c:(.text+0xa11): undefined reference to `gk_free'
/usr/bin/ld: graph.c:(.text+0xa29): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(graph.c.o): in function `libmetis__graph_ReadFromDisk':
graph.c:(.text+0xb11): undefined reference to `gk_rmpath'
/usr/bin/ld: graph.c:(.text+0xb7c): undefined reference to `gk_rmpath'
/usr/bin/ld: graph.c:(.text+0xb9c): undefined reference to `gk_errexit'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(graph.c.o): in function `libmetis__FreeSData':
graph.c:(.text+0x66a): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(initpart.c.o): in function `libmetis__Init2WayPartition':
initpart.c:(.text+0xb09): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: initpart.c:(.text+0xb42): undefined reference to `gk_errexit'
/usr/bin/ld: initpart.c:(.text+0xb59): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(initpart.c.o): in function `libmetis__GrowBisectionNode':
initpart.c:(.text+0xd1f): undefined reference to `gk_malloc'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(initpart.c.o): in function `libmetis__InitSeparator':
initpart.c:(.text+0x10fb): undefined reference to `gk_errexit'
/usr/bin/ld: initpart.c:(.text+0x1172): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: initpart.c:(.text+0x11b1): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(initpart.c.o): in function `GrowBisectionNode2':
initpart.c:(.text+0x12a1): undefined reference to `gk_malloc'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(kwayfm.c.o): in function `libmetis__Greedy_KWayOptimize':
kwayfm.c:(.text+0x71e1): undefined reference to `gk_errexit'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(sfm.c.o): in function `libmetis__FM_2WayNodeRefine1Sided':
sfm.c:(.text+0x190e): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: sfm.c:(.text+0x1942): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: sfm.c:(.text+0x198d): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: sfm.c:(.text+0x19b1): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: sfm.c:(.text+0x19f3): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(sfm.c.o):sfm.c:(.text+0x1a6b): more undefined references to `gk_CPUSeconds' follow
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(srefine.c.o): in function `libmetis__Allocate2WayNodePartitionMemory':
srefine.c:(.text+0x73): undefined reference to `gk_malloc'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(srefine.c.o): in function `libmetis__Refine2WayNode':
srefine.c:(.text+0x2ab): undefined reference to `gk_errexit'
/usr/bin/ld: srefine.c:(.text+0x314): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: srefine.c:(.text+0x361): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: srefine.c:(.text+0x380): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: srefine.c:(.text+0x3a1): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: srefine.c:(.text+0x3c1): undefined reference to `gk_CPUSeconds'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(srefine.c.o):srefine.c:(.text+0x3ff): more undefined references to `gk_CPUSeconds' follow
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(mincover.c.o): in function `libmetis__MinCover_Decompose':
mincover.c:(.text+0x403): undefined reference to `gk_free'
/usr/bin/ld: /home/z651w035/local/metis/lib/libmetis.a(mincover.c.o): in function `libmetis__MinCover':
mincover.c:(.text+0x8df): undefined reference to `gk_free'
collect2: error: ld returned 1 exit status
make[2]: *** [CMakeFiles/zhijian.dir/build.make:155: zhijian] Error 1
make[1]: *** [CMakeFiles/Makefile2:111: CMakeFiles/zhijian.dir/all] Error 2
make: *** [Makefile:136: all] Error 2
