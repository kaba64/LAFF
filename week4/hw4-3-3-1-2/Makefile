CC := gcc
LFLAGS = -lflame -lopenblas
TARGET_EXEC := out
OBJS = Symv_u_unb_var2.o axpy.o AXPY_Row_to_Column.o

$(TARGET_EXEC): $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)
%.o: %.c
	$(CC) -c $< -o $@
clean:
	rm $(OBJS) $(TARGET_EXEC) *~
run:
	@./$(TARGET_EXEC)
