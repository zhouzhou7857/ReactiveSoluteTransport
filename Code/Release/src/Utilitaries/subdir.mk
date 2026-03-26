################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Utilitaries/FluxPoint2D.cpp \
../src/Utilitaries/Laplace_Stehfest.cpp \
../src/Utilitaries/LinearSystem.cpp \
../src/Utilitaries/Point_Cgal.cpp \
../src/Utilitaries/RandomNumber.cpp \
../src/Utilitaries/Segment.cpp \
../src/Utilitaries/Structures.cpp \
../src/Utilitaries/UblasStructures.cpp \
../src/Utilitaries/scale.cpp 

OBJS += \
./src/Utilitaries/FluxPoint2D.o \
./src/Utilitaries/Laplace_Stehfest.o \
./src/Utilitaries/LinearSystem.o \
./src/Utilitaries/Point_Cgal.o \
./src/Utilitaries/RandomNumber.o \
./src/Utilitaries/Segment.o \
./src/Utilitaries/Structures.o \
./src/Utilitaries/UblasStructures.o \
./src/Utilitaries/scale.o 

CPP_DEPS += \
./src/Utilitaries/FluxPoint2D.d \
./src/Utilitaries/Laplace_Stehfest.d \
./src/Utilitaries/LinearSystem.d \
./src/Utilitaries/Point_Cgal.d \
./src/Utilitaries/RandomNumber.d \
./src/Utilitaries/Segment.d \
./src/Utilitaries/Structures.d \
./src/Utilitaries/UblasStructures.d \
./src/Utilitaries/scale.d 


# Each subdirectory must supply rules for building sources it contributes
src/Utilitaries/%.o: ../src/Utilitaries/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -I /home/zhouw/softs/ -I /home/zhouw/softs/boost/py37/gcc95/1.76.0/include/ -I /home/zhouw/softs/cgal/gcc95/5.5.2/include/ -I ./../../../Library -I /home/zhouw/softs/suitesparse/gcc95/5.10.1/include/ -I /usr/include/GL/ -frounding-math -O3 -DNDEBUG -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


