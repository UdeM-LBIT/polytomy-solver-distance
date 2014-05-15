g++ -shared -Wl,-soname,libpolytomysolver -o libpolytomysolver.so *.o -lpython2.7 -lboost_python -fPIC
