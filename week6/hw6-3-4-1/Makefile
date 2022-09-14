CC := gcc 
LFLAGS = -lflame -lopenblas
TARGET_EXEC := out
OBJS = solver.o axpy-gaussian.o axpy-gaussian-rank1.o ger.o dot-gaussian.o


$(TARGET_EXEC): $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)
%.o: %.c
	$(CC) -c $< -o $@
clean:
	rm $(OBJS) $(TARGET_EXEC) *~
run:
	@./$(TARGET_EXEC)
