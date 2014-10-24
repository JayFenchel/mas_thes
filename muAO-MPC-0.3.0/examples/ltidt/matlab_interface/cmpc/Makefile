CC = gcc 
FLAGS = -Os -Wall -Wstrict-prototypes
OPT = -O3 -funroll-loops
STD = -std=c99

all: libcmpc mpc_const.o mpc.o mpc_ref.o mpc_inc.o mpc_stc.o mtx_ops.o

libcmpc: mtx_ops.o mpc_const.o mpc.o mpc_ref.o mpc_stc.o mpc_inc.o
	ar rcs libcmpc.a mtx_ops.o mpc_const.o mpc.o mpc_ref.o mpc_inc.o mpc_stc.o

mpc_const.o: mpc_const.c
	$(CC) $(FLAGS) $(OPT) $(STD) -I./include -c mpc_const.c

mtx_ops.o: mtx_ops.c
	$(CC) $(FLAGS) $(OPT) $(STD) -I./include -c mtx_ops.c

mpc.o: mpc.c
	$(CC) $(FLAGS) $(OPT) $(STD) -I./include -c mpc.c

mpc_ref.o: mpc_ref.c
	$(CC) $(FLAGS) $(OPT) $(STD) -I./include -c mpc_ref.c

mpc_inc.o: mpc_inc.c
	$(CC) $(FLAGS) $(OPT) $(STD) -I./include -c mpc_inc.c

mpc_stc.o: mpc_stc.c
	$(CC) $(FLAGS) $(OPT) $(STD) -I./include -c mpc_stc.c

clean:
	rm *.o libcmpc.a
