CC := gcc
LFLAGS = -lflame -lopenblas
TARGET_EXEC := out
SRCS := Mvmult_t_n_var2.c
OBJS = Mvmult_t_n_var2.o AXPY_Row_to_Column.o


$(TARGET_EXEC): $(OBJS)
	$(CC) $^ $(LFLAGS) -o $@
Mvmult_t_n_var2.o: Mvmult_t_n_var2.c
	$(CC) $(LFLAGS) -c $<
AXPY_Row_to_Column.o: AXPY_Row_to_Column.c
	$(CC) $(LFLAGS) -c $<
clean:
	rm $(OBJS) $(TARGET_EXEC) *~
run:
	@./$(TARGET_EXEC)
