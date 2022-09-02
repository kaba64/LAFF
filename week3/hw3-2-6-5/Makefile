CC := gcc
LFLAGS = -lflame -lopenblas
TARGET_EXEC := out
SRCS := Transpose.c
OBJS=$(SRCS:.c=.o)


$(TARGET_EXEC): $(OBJS)
	$(CC) $^ $(LFLAGS) -o $@
$(OBJS): $(SRCS)
	$(CC) $(LFLAGS) -c $<
clean:
	rm $(OBJS) $(TARGET_EXEC) *~
run:
	@./$(TARGET_EXEC)
