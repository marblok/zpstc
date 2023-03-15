# ZPSTC

**DSPElib applications used during ZPSTC laboratory**

The files are prepared for the use with Visusal Studio Code, MSys2 with MinGW x64 and gcc compiler from MSys2.

1. The following directories structure is assumed:

DSPElib library folder (needed only for library compilation with task „Install DSPElib”, https://git.pg.edu.pl/dspe)
	MAIN_DIR\DSPElib
Folder with compiled DSPElib prepared for linking (e.g. with subfolder x86_64-w64-mingw32-gcc_12.2.0)
	MAIN_DIR\_DSPE_lib_minGW_
ZPSTC main folder
	MAIN_DIR\ZPSTC
ZPSTC excercises folders
	MAIN_DIR\ZPSTC\Ex1
	MAIN_DIR\ZPSTC\Ex2
	...
	
2. Configuration files in .vscode subfolder that might need adjusting (per excercise) depending on you environment setup:
 - c_cpp_properties.json (utilized by Intellisense)
   - DSPElib subfolder in "includePath" per configuration
   ~~~
	"${workspaceFolder}/../../_DSPE_lib_minGW_/x86_64-w64-mingw32-gcc_12.2.0/include/**"
   ~~~
 - tasks.json (compilation tasks: menu: Terminal > Run Tasks...)
	- "env" in "windows" sections: path to msys2 tools (per task)
 - launch.json
    - "miDebuggerPath": path to gdb if it's not in system path
	
