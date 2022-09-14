CC := gcc 
LFLAGS = -lflame -lopenblas
TARGET_EXEC := out
OBJS = Utrsv_unb_var1.o dot-gaussian.o

$(TARGET_EXEC): $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)
%.o: %.c
	$(CC) -c $< -o $@
clean:
	rm $(OBJS) $(TARGET_EXEC) *~
run:
	@./$(TARGET_EXEC)
