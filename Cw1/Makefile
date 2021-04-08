# Run: make Release 
# Run: make Debug 
CC=g++
# comflag = -m32
# comflag = -m64
comflag = $(COMFLAG)

DSPElib_DIR = ../../_DSPE_lib_minGW_/$(DSPElib_SUBDIR)

SRC_DIR = .
SRC_CPP_SUBDIR = .

#DFLAGS         = -DWIN32 -DDEVCPP 

# -D INCLUDE_DSPE_EXAMPLES # TODO: użycie w ramach kompilacji Main.cpp w trybie DEBUG
ifeq ($(MODE),Release)
	CFLAGS = $(comflag) -std=c++0x -O3 -Wall -c -fmessage-length=0 -fno-strict-aliasing 
	LINKER_FLAGS = $(comflag)  -s -static-libgcc -static-libstdc++ $(MISC_LINKER_FLAGS)
	INCLUDES := -I"$(DSPElib_DIR)/include" -I"$(DSPElib_DIR)/include/rls"
	DSPElib_FULLDIR = $(DSPElib_DIR)/rls
	EXE_FILENAME = Cw1_zad3_rls.exe
else
	CFLAGS   = $(comflag) -std=c++0x -O0 -g3 -Wall -c -fmessage-length=0 -W -Wshadow -Wconversion -fstrict-aliasing -fmax-errors=5
	LINKER_FLAGS   = $(comflag)  -static-libgcc -static-libstdc++ $(MISC_LINKER_FLAGS)
	INCLUDES := -I"$(DSPElib_DIR)/include" -I"$(DSPElib_DIR)/include/dbg"
	DSPElib_FULLDIR = $(DSPElib_DIR)/dbg
	EXE_FILENAME = Cw1_zad3_dbg.exe
endif
# -U__STRICT_ANSI__ jest potrzebne do kompilacji debug_new.cpp, jezeli pominac ten plik to mozna rowniez wyrzucic te opcje
#CFLAGS_debug   = $(comflag) -std=c++0x -O0 -g3 -Wall -c -fmessage-length=0 -W -Wshadow -Wco#nversion -fstrict-aliasing -U__STRICT_ANSI__

SOURCES_NAMES = 
SOURCES_NAMES += Cw1_zad3.cpp 
SOURCES = $(addprefix $(SRC_CPP_SUBDIR)/,$(SOURCES_NAMES))

SOURCES_DBG =
# SOURCES_DBG += $(SRC_DIR)/Main.cpp

# ################################################# #
# DEBUG
OBJECTS 		  := $(SOURCES:%.cpp=$(OUT_DIR)/%.o)
DEPENDS 		  := $(SOURCES:%.cpp=$(OUT_DIR)/%.d)

# ################################################# #
-include $(DEPENDS)

all: build


# ########################################################################################### #	
# ########################################################################################### #	
build: $(SRC_DIR)/$(EXE_FILENAME)

$(SRC_DIR)/$(EXE_FILENAME): $(OBJECTS)
	@echo $(EXE_FILENAME)
	$(CC) -L$(DSPElib_FULLDIR) $(OBJECTS) -o"$(SRC_DIR)/$(EXE_FILENAME)" $(LINKER_FLAGS) -lDSPE $(LIBS)

# ########################################################################################### #	
# ########################################################################################### #	
# Z podanej listy usuwany $(OUT_DIR_WIN_RLS)/ oraz '.o' zamieniamy na '.cpp'
$(OBJECTS): $(OUT_DIR)/%.o : %.cpp
	@echo $(@D) $< $@

	#mkdir -p $(OUT_DIR)/$(SRC_CPP_SUBDIR)
	mkdir -p $(@D)
	$(CC) $(DFLAGS) $(CFLAGS) $(INCLUDES) -MMD $< -o $@


clean:
	@echo MODE: $(MODE)

	@if [ -d "$(OUT_DIR)" ]; then \
		echo "cleaning $(OUT_DIR_DBG) ..."; \
		#find $(OUT_DIR)/ -name "*.o" -type f -delete; \
		rm -rf $(OUT_DIR)/*.d; \
		rm -rf $(OUT_DIR)/*.o; \
		rm -rf $(OUT_DIR); \
		echo "cleaned $(OUT_DIR)"; \
	fi 
	rm -rf "$(SRC_DIR)/$(EXE_FILENAME)"; \
	#rm -rf "$(SRC_DIR)/*.gif"; \
	#rm -rf "$(SRC_DIR)/*.dot"; \
	
.PHONY: all build clean

