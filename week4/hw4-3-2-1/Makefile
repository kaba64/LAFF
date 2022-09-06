CC := gcc
LFLAGS = -lflame -lopenblas
TARGET_EXEC := out
SRCS := Trmvp_un_unb_var1.c
OBJS = Trmvp_un_unb_var1.o dot.o


$(TARGET_EXEC): $(OBJS)
	$(CC) $^ $(LFLAGS) -o $@
Trmvp_un_unb_var1.o: Trmvp_un_unb_var1.c
	$(CC) $(LFLAGS) -c $<
dot.o: dot.c
	$(CC) $(LFLAGS) -c $<
clean:
	rm $(OBJS) $(TARGET_EXEC) *~
run:
	@./$(TARGET_EXEC)
