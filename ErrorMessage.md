[ 84%] Building CXX object CMakeFiles/zhijian.dir/src/main.cpp.o
/home/z651w035/codes/Zhijian/src/main.cpp: In function ‘zhijian::SimParams parseConfig(const string&)’:
/home/z651w035/codes/Zhijian/src/main.cpp:78:36: error: variable ‘std::istringstream iss’ has initializer but incomplete type
   78 |         std::istringstream iss(line);
      |                                    ^
/home/z651w035/codes/Zhijian/src/main.cpp: In function ‘int main(int, char**)’:
/home/z651w035/codes/Zhijian/src/main.cpp:307:31: error: ‘setw’ is not a member of ‘std’; did you mean ‘set’?
  307 |             std::cout << std::setw(6) << iter << "  "
      |                               ^~~~
      |                               set
/home/z651w035/codes/Zhijian/src/main.cpp:308:50: error: ‘setprecision’ is not a member of ‘std’
  308 |                       << std::scientific << std::setprecision(4)
      |                                                  ^~~~~~~~~~~~
/home/z651w035/codes/Zhijian/src/main.cpp:309:31: error: ‘setw’ is not a member of ‘std’; did you mean ‘set’?
  309 |                       << std::setw(12) << time << "  "
      |                               ^~~~
      |                               set
/home/z651w035/codes/Zhijian/src/main.cpp:310:31: error: ‘setw’ is not a member of ‘std’; did you mean ‘set’?
  310 |                       << std::setw(12) << residual << "\n";
      |                               ^~~~
      |                               set
/home/z651w035/codes/Zhijian/src/main.cpp:320:32: error: aggregate ‘std::ostringstream fname’ has incomplete type and cannot be defined
  320 |             std::ostringstream fname;
      |                                ^~~~~
/home/z651w035/codes/Zhijian/src/main.cpp:322:34: error: ‘setfill’ is not a member of ‘std’; did you mean ‘fill’?
  322 |                   << "_" << std::setfill('0') << std::setw(4) << output_count << ".vtu";
      |                                  ^~~~~~~
      |                                  fill
/home/z651w035/codes/Zhijian/src/main.cpp:322:55: error: ‘setw’ is not a member of ‘std’; did you mean ‘set’?
  322 |                   << "_" << std::setfill('0') << std::setw(4) << output_count << ".vtu";
      |                                                       ^~~~
      |                                                       set
make[2]: *** [CMakeFiles/zhijian.dir/build.make:76: CMakeFiles/zhijian.dir/src/main.cpp.o] Error 1
make[1]: *** [CMakeFiles/Makefile2:111: CMakeFiles/zhijian.dir/all] Error 2
make: *** [Makefile:136: all] Error 2

