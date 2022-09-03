CC := gcc
LFLAGS = -lflame -lopenblas
TARGET_EXEC := out
SRCS := Mvmult_n_var2.c
OBJS = Mvmult_n_var2.o axpy.o


$(TARGET_EXEC): $(OBJS)
	$(CC) $^ $(LFLAGS) -o $@
Mvmult_n_var2.o: Mvmult_n_var2.c
	$(CC) $(LFLAGS) -c $<
axpy.o: axpy.c
	$(CC) $(LFLAGS) -c $<
clean:
	rm $(OBJS) $(TARGET_EXEC) *~
run:
	@./$(TARGET_EXEC)
