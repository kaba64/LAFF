CC := gcc
LFLAGS = -lflame -lopenblas
TARGET_EXEC := out
SRCS := Mvmult_n_var1.c
OBJS = Mvmult_n_var1.o dot.o


$(TARGET_EXEC): $(OBJS)
	$(CC) $^ $(LFLAGS) -o $@
Mvmult_n_var1.o: Mvmult_n_var1.c
	$(CC) $(LFLAGS) -c $<
dot.o: dot.c
	$(CC) $(LFLAGS) -c $<
clean:
	rm $(OBJS) $(TARGET_EXEC) *~
run:
	@./$(TARGET_EXEC)
