CC := gcc
LFLAGS = -lflame -lopenblas
TARGET_EXEC := out
SRCS := Mvmult_axpy.c
OBJS = Mvmult_axpy.o axpy.o


$(TARGET_EXEC): $(OBJS)
	$(CC) $^ $(LFLAGS) -o $@
Mvmult_axpy.o: Mvmult_axpy.c
	$(CC) $(LFLAGS) -c $<
axpy.o: axpy.c
	$(CC) $(LFLAGS) -c $<
clean:
	rm $(OBJS) $(TARGET_EXEC) *~
run:
	@./$(TARGET_EXEC)
