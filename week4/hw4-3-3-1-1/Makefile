CC := gcc
LFLAGS = -lflame -lopenblas
TARGET_EXEC := out
OBJS = Symv_u_unb_var1.o dot.o Dot_Both_column_vector.o

$(TARGET_EXEC): $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)
%.o: %.c
	$(CC) -c $< -o $@
clean:
	rm $(OBJS) $(TARGET_EXEC) *~
run:
	@./$(TARGET_EXEC)
