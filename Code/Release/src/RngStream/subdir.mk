################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/RngStream/rngstream.cpp 

OBJS += \
./src/RngStream/rngstream.o 

CPP_DEPS += \
./src/RngStream/rngstream.d 


# Each subdirectory must supply rules for building sources it contributes
src/RngStream/%.o: ../src/RngStream/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I ./../Library -I /usr/include/suitesparse/ -I /usr/include/GL/ -frounding-math -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


