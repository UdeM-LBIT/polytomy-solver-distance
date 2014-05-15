The source is given as a Qt project, though there are no dependencies to any of the Qt libraries.

However, whichever environment is used, the user must ensure that the following Qt
compiler directives are handled : 
DEFINES -= UNICODE
DEFINES += _MBCS
QMAKE_CXXFLAGS -= DUNICODE
QMAKE_CXXFLAGS += -std=c++0x

Essentially, the project uses multibyte string encoding, and requires the use of c++0x (for the unordered_set class).

If Qt is installed, the simplest to to get into the folder, and enter the two commands
> qmake PolytomySolver.pro "CONFIG+=RELEASE"
> make

Otherwise, the provided Makefile might be enough to simply run 
> make

And otherwise, the simplest way is to include all .h and .cpp in an Eclipse or QCreator project (or even VC++), and build.