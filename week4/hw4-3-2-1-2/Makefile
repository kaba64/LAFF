CC := gcc
LFLAGS = -lflame -lopenblas
TARGET_EXEC := out
SRCS := Trmvp_un_unb_var2.c
OBJS = Trmvp_un_unb_var2.o axpy.o


$(TARGET_EXEC): $(OBJS)
	$(CC) $^ $(LFLAGS) -o $@
Trmvp_un_unb_var2.o: Trmvp_un_unb_var2.c
	$(CC) $(LFLAGS) -c $<
axpy.o: axpy.c
	$(CC) $(LFLAGS) -c $<
clean:
	rm $(OBJS) $(TARGET_EXEC) *~
run:
	@./$(TARGET_EXEC)
