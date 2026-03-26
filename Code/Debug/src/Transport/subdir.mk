################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Transport/Particle.cpp \
../src/Transport/Projection.cpp \
../src/Transport/Transfer.cpp \
../src/Transport/Transport.cpp 

OBJS += \
./src/Transport/Particle.o \
./src/Transport/Projection.o \
./src/Transport/Transfer.o \
./src/Transport/Transport.o 

CPP_DEPS += \
./src/Transport/Particle.d \
./src/Transport/Projection.d \
./src/Transport/Transfer.d \
./src/Transport/Transport.d 


# Each subdirectory must supply rules for building sources it contributes
src/Transport/%.o: ../src/Transport/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I/usr/include -I/usr/include/c++/5/bits -I/home/marwa/libraries/boost_1_61_0 -I/home/marwa/libraries/boost-numeric-bindings -I/home/marwa/libraries/SuiteSparse/UMFPACK/Include -I/home/marwa/libraries/SuiteSparse/AMD/Include -I/home/marwa/libraries/CGAL-4.12-beta1/include -I/home/marwa/libraries/SuiteSparse/include -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


