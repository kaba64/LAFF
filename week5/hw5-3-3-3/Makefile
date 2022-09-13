CC := gcc 
LFLAGS = -lflame -lopenblas
TARGET_EXEC := out
OBJS = Gemm_unb_var2.o gevm_dot.o dot.o

$(TARGET_EXEC): $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)
%.o: %.c
	$(CC) -c $< -o $@
clean:
	rm $(OBJS) $(TARGET_EXEC) *~
run:
	@./$(TARGET_EXEC)
