1. Download the project with `git clone --recursive https://github.com/Vekleon/sand-simulation.git`
2. Inside the project folder create a build folder:
			mkdir build
			cd build
			cmake ..
   If you are on Mac or Linux, then do:
			make
3. Run `cmake ..` in this directory. If on Windows, be sure to use `x86`. This should create a `.sln` file which can be opened with Visual Studio.
4. Running cmake should have also generated a `.exe` in either the Debug or Release folder, depending on the flags used with `cmake`. Move this to the `build` folder so that the file links are not broken.
5. Run the executable to see a bunny. :)