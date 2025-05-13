CC=riscv64-linux-gnu-gcc
AR=riscv64-linux-gnu-ar
ARCH=rv64gcv_zba

libnttvec.a: nttvec.o
	$(AR) rcs $@ $^

nttvec.o: src/ntt.s
	$(CC) -c -march=$(ARCH) $^ -o $@
