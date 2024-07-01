# Printing
println("Hello world!")

# Assigning variables
# Types do not need to be specified
myInt = 12
println(typeof(myInt))

#=
Multi-line commenting
=#

# Basic math
sum1 = 3 + 7
diff = 3 - 7
product = 3 * 7
quotient = 21 / 3
power = 3 ^ 3
mod = 23 % 7
println(power)

# Converting between types
myFloat = convert(Float64, myInt)
println(myFloat)

# Creating strings
myString = "Hello!"
myString2 = """This allows "quotation marks" within the string."""
println(myString2)

# String interpolation
a = 3
b = 7
math = "a = $a. b = $b. a + b = $(a+b)."
println(math)

# String concatenation
firstHalf = "Combine this half "
secondHalf = "and this half."
combined = string(firstHalf, secondHalf)
combined2 = firstHalf * secondHalf
combined3 = "$firstHalf$secondHalf"
println(combined3)

# Dictionaries are unordered sets of key-value pairs
capitals = Dict("USA" => "Washington, D.C.", "Canada" => "Ottawa")
capitals["Mexico"] = "Mexico City"
println(capitals["USA"])
println(pop!(capitals, "Canada"))

# Tuples are immutable, ordered lists
# Indexing starts from 1
colors = ("red", "yellow", "blue")
println(colors[2])

# Arrays are mutable, ordered lists
fibonacci = [1, 1, 2, 3, 5, 8]
push!(fibonacci, 13)
println(pop!(fibonacci))
twod = [[1, 3, 5], [2, 4, 6]]
println(twod[1][2])

# while loops
n = 1
while n < 1000
    global n *= 2
    println(n)
end

# for loops
for n in 1:10
    println("$(2^n)")
end

addition = zeros(5, 5)
for i in 1:5
    for j in 1:5
        addition[i, j] = i + j
    end
end
println(addition[3, 5])

multiplication = zeros(5,5)
for i in 1:5, j in 1:5
    multiplication[i, j] = i * j
end
println(multiplication[3, 5])

powers = [i ^ j for i in 1:5, j in 1:3]
println(powers[4, 2])

# Conditionals
grade = 40
if grade >= 70
    println("You passed!")
elseif grade >= 60
    println("You got a D.")
else
    println("You failed.")
end

(a > b) ? println(a) : println(b)

# Functions
function square(x)
    x ^ 2
end
println(square(3))

cube(x) = x ^ 3
println(cube(3))

sqrt = x -> x ^ 0.5
println(sqrt(25))

# Mutating vs non-mutating functions
# Mutating functions are followed by an !
unordered = [2, 3, 1]
println(sort(unordered))
println(unordered)
println(sort!(unordered))
println(unordered)

# Broadcasting functions
# Broadcasting functions are followed by a .
# Broadcasting functions apply a function to each element of an input 
# rather than to the input as a whole.
println(square(addition))
println(square.(addition))

# Packages
# using Pkg
# Pkg.add("Plots")
using Plots

# Plotting
x = 0:0.1:(2*π)
y = sin.(x)
gr()
plot(x, y, label = "lines")
scatter!(x, y, label = "points")
display(title!("Plot of sin(x)"))

p1 = plot(x, sin.(x))
p2 = plot(x, cos.(x))
plot(p1, p2, layout = (1,2), legend = false)

# Multiple dispatch
# One function can have different purposes
# depending on the types of the input variables.
println(3 + 7)
println(3 + 7.0)
import Base: +
+(x::String, y::String) = string(x, y)
println(firstHalf + secondHalf)

# Structs
# Similar to classes in Java
struct Person
    firstName
    lastName
end
me = Person("Andrew","Vodinh-Ho")
println(me.firstName)

mutable struct Profile
    name::String
    age::Float64
    isAlive::Bool

    function Profile(name, age)
        new(name, age, true)
    end
end
me2 = Profile("Andrew", 19)

function birthday(person::Profile)
    person.age+=1
end
birthday(me2)
println(me2.age)

# Testing speed
#randoms = rand(10^7)
#Pkg.add("BenchmarkTools")
#using BenchmarkTools
#bench = @benchmark sum($randoms)
#println(minimum(bench.times))

# Linear algebra
matrixA = rand(1:4, 3, 3)
matrixB = matrixA
matrixC = copy(matrixA)
matrixA[1, 1]=17
println(matrixB)
println(matrixC)
conjugateTranspose = matrixA'
rhs = ones(3)
sol = matrixA\rhs
println(sol)

