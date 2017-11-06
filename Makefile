
# This is the makefile for APA
# 1. Use g++ to compile
# 2. Link libraray Eign to solve linear equations
# 3. Compile under c++11 standard
# 4. After compling, all the *.o files will be located in ./obj
# 5. After compling, the executable file be located in ./

TARGET = APA	

.PHONY:clean

CC = g++
FLAGS = -I ./Eigen -stdlib=libc++ -std=c++11
CPP_FILES = $(shell ls *.cpp)
BASE = $(basename $(CPP_FILES))
OBJDIR = obj
OBJS = $(addsuffix .o, $(addprefix $(OBJDIR)/,$(BASE)))

all: $(TARGET)
	@echo "Finish compiling..."
$(OBJS): | $(OBJDIR)
$(OBJDIR):
	@mkdir -p $(OBJDIR)

$(TARGET):$(OBJS) $(OBJDIR)
	@${CC} ${FLAGS} -o $(TARGET) $(OBJS)

$(OBJDIR)/%.o:%.cpp
	@${CC} ${FLAGS} -c -o $@ $<
	@echo "Compiling $<......"
clean:
	@rm -r $(OBJDIR) $(TARGET)