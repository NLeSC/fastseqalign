CC = gcc
CXX = g++
CFLAGS := -g -Wall -pipe -O2 -mavx2 #-g -fsanitize=address
CXXFLAGS := $(CFLAGS)
LOBJS = ssw.o
LCPPOBJS = ssw_cpp.o
PROG = ssw_test
LIB = libssw.so
EXAMPLE = example_c
EXAMPLE_CPP = example_cpp
JAVA_JAR = ssw.jar
JAVA_LIB = libsswjni.so
C_SRC_PATH = ../ssw_native
JAVA_CLASS_PATH = ./nl/escience/alignment
JAVA_INLCUDES = -I"$(JAVA_HOME)/include" -I"$(JAVA_HOME)/include/linux" -I/usr/lib/jvm/java-11-openjdk-amd64/include -I/usr/lib/jvm/java-1.11.0-openjdk-amd64/include/linux
JAVA_OBJ =  $(JAVA_CLASS_PATH)/Alignment.class $(JAVA_CLASS_PATH)/LocalSequenceAlignment.class $(JAVA_CLASS_PATH)/BoundedLocalSequenceAlignment.class   $(JAVA_CLASS_PATH)/Example.class

.PHONY: all default java clean

core: $(PROG)

default: $(PROG) $(EXAMPLE) $(EXAMPLE_CPP) $(LIB) 

all: default java

java: $(JAVA_JAR) $(JAVA_LIB)

$(LIB): ssw.c ssw.h
	$(CC) $(CFLAGS) -fPIC -shared -rdynamic -o $@ $<

$(PROG): main.c kseq.h

$(EXAMPLE): example.c

ifdef __arm__ # (M1)
$(PROG) $(EXAMPLE): $(LOBJS)
	$(CC) -o $@ $(filter-out %.h,$^) $(CFLAGS) -lm -lz -march=armv8-a+fp+simd+crypto+crc
else # x86(Intel)
$(PROG) $(EXAMPLE): $(LOBJS)
	$(CC) -o $@ $(filter-out %.h,$^) $(CFLAGS) -lm -lz
endif


$(EXAMPLE_CPP): example.cpp $(LOBJS) $(LCPPOBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS) -lm -lz

$(JAVA_LIB): $(C_SRC_PATH)/sswjni.c $(C_SRC_PATH)/ssw.c $(C_SRC_PATH)/ssw.h
	$(CC) $(CFLAGS) $(JAVA_INLCUDES) -fPIC -shared -rdynamic -o $@ $< $(C_SRC_PATH)/ssw.c

$(JAVA_JAR): $(JAVA_OBJ)
	jar cvfe $@ nl.escience.alignment.Example $^

%.class: %.java
	javac -cp  ./ -Xlint $<
	
ssw.o: ssw.c ssw.h
	$(CC) -c -o $@ $< $(CFLAGS)

ssw_cpp.o: ssw_cpp.cpp ssw_cpp.h ssw.h
	$(CXX) -c -o $@ $< $(CXXFLAGS)

clean:
	-rm -f $(LOBJS) $(LCPPOBJS) $(PROG) $(LIB) $(EXAMPLE) $(EXAMPLE_CPP) $(JAVA_LIB) $(JAVA_JAR) $(JAVA_OBJ) *~ 
