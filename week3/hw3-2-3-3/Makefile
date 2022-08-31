CC := gcc
LFLAGS = -lflame -llapack
TARGET_EXEC := out
SRCS := Set_to_diagonal_matrix.c
OBJS=$(SRCS:.c=.o)


$(TARGET_EXEC): $(OBJS)
	$(CC) $^ $(LFLAGS) -o $@
$(OBJS): $(SRCS)
	$(CC) $(LFLAGS) -c $<
clean:
	rm $(OBJS) $(TARGET_EXEC) *~
run:
	@./$(TARGET_EXEC)
