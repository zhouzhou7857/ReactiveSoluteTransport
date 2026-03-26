################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Visualisation/DFNVisu.cpp \
../src/Visualisation/DisplayResults.cpp 

OBJS += \
./src/Visualisation/DFNVisu.o \
./src/Visualisation/DisplayResults.o 

CPP_DEPS += \
./src/Visualisation/DFNVisu.d \
./src/Visualisation/DisplayResults.d 


# Each subdirectory must supply rules for building sources it contributes
src/Visualisation/%.o: ../src/Visualisation/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/include -I/usr/include/c++/5/bits -I/home/marwa/libraries/boost_1_61_0 -I/home/marwa/libraries/boost-numeric-bindings -I/home/marwa/libraries/SuiteSparse/UMFPACK/Include -I/home/marwa/libraries/SuiteSparse/AMD/Include -I/home/marwa/libraries/CGAL-4.12-beta1/include -I/home/marwa/libraries/SuiteSparse/include -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


