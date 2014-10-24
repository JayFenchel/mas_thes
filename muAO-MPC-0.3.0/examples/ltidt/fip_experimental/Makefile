CC = gcc
FLAGS = -Os -Wall -Wstrict-prototypes -pedantic
OPT = -O3 -funroll-loops
STD = -std=c89

all: lcmpc motor

lcmpc:
	make -C cmpc

motor: main_motor.c lcmpc
	$(CC) $(FLAGS) $(OPT) $(STD) main_motor.c -Lcmpc -lcmpc -o motor

clean:
	rm -rf cmpc python __pycache__ motor
