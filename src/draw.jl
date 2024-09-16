struct SquarePointLength
    p::Point
    l::Real
end

SquarePointLength(x::Int, y::Int, l::Real) = SquarePointLength(Point(x,y), l)
SquarePointLength(p::CartesianIndex{2}, l::Real) = SquarePointLength(Point(p), l)

function draw!(img::AbstractArray{T, 2}, square::SquarePointLength, color::T) where {T<:Colorant}
    ind = axes(img)
    x = square.p.x; y = square.p.y
    l = square.l
    (first(ind[2]) <= x <= last(ind[2]) - l && first(ind[1]) <= y <= last(ind[1]) - l) || error("Square is out of the bounds of image")
    draw!(img, Rectangle(Point(x, y), l, l), color)
end

function draw!(img::AbstractArray{T, 2}, square::SquarePointLength, color::T) where {T<:Colorant}
    ind = axes(img)
    x = square.p[1]; y = square.p[2]
    l = square.l
    (first(ind[2]) <= x <= last(ind[2]) - l && first(ind[1]) <= y <= last(ind[1]) - l) || error("Square is out of the bounds of image")
    draw!(img, Rectangle(Point(x, y), l, l), color)
end

# Function to draw a square
function draw_square!(image, x, y, size)
    image[x:x+size, y:y+size] .= RGB(1, 0, 0)  # Red square
end

# Function to draw a triangle
function draw_triangle!(image, x, y, size)
    for i in x:x+size
        for j in y:y+size
            if i - x + j - y <= size
                image[i, j] = RGB(0, 0, 1)  # Blue triangle
            end
        end
    end
end

function draw_diamond!(image, x, y, size)
    for i in x:x+size
        for j in y:y+size
            if abs(i - x - size//2) + abs(j - y - size//2) <= size//2
                image[i, j] = RGB(0, 1, 0)  # Green diamond
            end
        end
    end
end

