CC := gcc 
LFLAGS = -lflame -lopenblas
TARGET_EXEC := out
OBJS = GJ_Inverse_alt_unb.o ger.o ax.o CV_Set_zero_unb.o Rank1-update.o Scale.o RV_Scale_unb.o

$(TARGET_EXEC): $(OBJS)
	$(CC) -o $@ $^ $(LFLAGS)
%.o: %.c
	$(CC) -c $< -o $@
clean:
	rm $(OBJS) $(TARGET_EXEC) *~
run:
	@./$(TARGET_EXEC)
