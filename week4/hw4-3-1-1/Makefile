CC := gcc
LFLAGS = -lflame -lopenblas
TARGET_EXEC := out
SRCS := Mvmult_t_n_var1.c
OBJS = Mvmult_t_n_var1.o Dot_Both_column_vector.o


$(TARGET_EXEC): $(OBJS)
	$(CC) $^ $(LFLAGS) -o $@
Mvmult_t_n_var1.o: Mvmult_t_n_var1.c
	$(CC) $(LFLAGS) -c $<
Dot_Both_column_vector.o: Dot_Both_column_vector.c
	$(CC) $(LFLAGS) -c $<
clean:
	rm $(OBJS) $(TARGET_EXEC) *~
run:
	@./$(TARGET_EXEC)
