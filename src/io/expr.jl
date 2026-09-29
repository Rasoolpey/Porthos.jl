# Safe reader for parameter expressions such as "2.0 * M_PI * 60.0".
#
# PHPS pastes these strings into generated C++ (and occasionally Python), so they follow C
# arithmetic: + - * / with the usual precedence, left associative, unary signs, parentheses,
# decimal literals and the constant M_PI. The reader evaluates in Float64 in the same order
# as C, so the result is bit-identical to what PHPS computed. Nothing is ever passed to
# `eval`; anything outside this grammar is an error.
#
#   expr   := term  (('+' | '-') term)*
#   term   := unary (('*' | '/') unary)*
#   unary  := ('+' | '-') unary | atom
#   atom   := number | constant | func '(' expr ')' | '(' expr ')'

"""
    ParamExprError

Raised when a parameter string is not a valid arithmetic expression.
"""
struct ParamExprError <: Exception
    expr::String
    msg::String
end
Base.showerror(io::IO, e::ParamExprError) =
    print(io, "ParamExprError: ", e.msg, " in parameter expression \"", e.expr, "\"")

const EXPR_CONSTANTS = Dict{String,Float64}(
    "M_PI" => Float64(π),
    "pi" => Float64(π),
    "M_E" => Float64(ℯ),
)

const EXPR_FUNCTIONS = Dict{String,Function}(
    "sqrt" => sqrt,
)

mutable struct _ExprLexer{T}
    s::String
    pos::Int                           # byte index of the next unread character
    env::Dict{String,T}                # named values (parameters), checked before constants
end

function _skipspace!(lx::_ExprLexer)
    while lx.pos <= ncodeunits(lx.s) && isspace(lx.s[lx.pos])
        lx.pos = nextind(lx.s, lx.pos)
    end
end

_peek(lx::_ExprLexer) = (_skipspace!(lx); lx.pos <= ncodeunits(lx.s) ? lx.s[lx.pos] : '\0')

function _take!(lx::_ExprLexer, c::Char)
    _peek(lx) == c || throw(ParamExprError(lx.s, "expected '$c' at position $(lx.pos)"))
    lx.pos = nextind(lx.s, lx.pos)
end

const _NUMBER_RE = r"\G(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?"
const _IDENT_RE = r"\G[A-Za-z_][A-Za-z_0-9]*"

function _expr(lx::_ExprLexer)
    v = _term(lx)
    while true
        c = _peek(lx)
        if c == '+'
            _take!(lx, c); v = v + _term(lx)
        elseif c == '-'
            _take!(lx, c); v = v - _term(lx)
        else
            return v
        end
    end
end

function _term(lx::_ExprLexer)
    v = _unary(lx)
    while true
        c = _peek(lx)
        if c == '*'
            _take!(lx, c); v = v * _unary(lx)
        elseif c == '/'
            _take!(lx, c); v = v / _unary(lx)
        else
            return v
        end
    end
end

function _unary(lx::_ExprLexer)
    c = _peek(lx)
    c == '-' && (_take!(lx, c); return -_unary(lx))
    c == '+' && (_take!(lx, c); return _unary(lx))
    return _atom(lx)
end

function _atom(lx::_ExprLexer)
    c = _peek(lx)
    if c == '('
        _take!(lx, '(')
        v = _expr(lx)
        _take!(lx, ')')
        return v
    end
    m = match(_NUMBER_RE, lx.s, lx.pos)
    if m !== nothing
        lx.pos += ncodeunits(m.match)
        return _expr_number(lx, parse(Float64, m.match))
    end
    m = match(_IDENT_RE, lx.s, lx.pos)
    if m !== nothing
        name = String(m.match)
        lx.pos += ncodeunits(name)
        if haskey(lx.env, name)
            return lx.env[name]
        elseif haskey(EXPR_CONSTANTS, name)
            return _expr_number(lx, EXPR_CONSTANTS[name])
        elseif haskey(EXPR_FUNCTIONS, name)
            _take!(lx, '(')
            v = _expr(lx)
            _take!(lx, ')')
            return _expr_call(lx, EXPR_FUNCTIONS[name], v)
        end
        throw(ParamExprError(lx.s, "unknown name '$name'"))
    end
    c == '\0' && throw(ParamExprError(lx.s, "unexpected end of expression"))
    throw(ParamExprError(lx.s, "unexpected character '$c' at position $(lx.pos)"))
end

"""
    parse_param_expr(s::AbstractString, env = Dict{String,Float64}()) -> Float64

Evaluate a parameter expression (numbers, `+ - * /`, parentheses, `M_PI`, `pi`, `M_E`,
`sqrt`, and the names in `env`) with C evaluation order in `Float64`. Never calls `eval`.

```jldoctest
julia> Porthos.parse_param_expr("2.0 * M_PI * 60.0") == 2.0 * π * 60.0
true
```
"""
parse_param_expr(s::AbstractString, env::AbstractDict = Dict{String,Float64}()) =
    eval_param_expr(Float64, s, env)

"""
    eval_param_expr(T, s, env) -> T

`parse_param_expr` in the number type `T`: every literal, constant and `env` value is
converted to `T` (a literal is first read as the nearest Float64, as the model reads its
parameters), and the arithmetic is `T`'s. With an interval type this gives an
outward-rounded enclosure of the Float64 expression's exact value (the ROA containment
audit).
"""
function eval_param_expr(::Type{T}, s::AbstractString, env::AbstractDict) where {T}
    lx = _ExprLexer{T}(String(s), 1, Dict{String,T}(string(k) => T(v) for (k, v) in env))
    _peek(lx) == '\0' && throw(ParamExprError(lx.s, "empty expression"))
    v = _expr(lx)
    _peek(lx) == '\0' || throw(ParamExprError(lx.s, "trailing input at position $(lx.pos)"))
    return v
end

_expr_number(::_ExprLexer{T}, v::Float64) where {T} = T(v)
_expr_call(::_ExprLexer{Float64}, f, v) = Float64(f(v))
_expr_call(::_ExprLexer, f, v) = f(v)

"""
    param_value(x) -> Float64

A numeric parameter as stored in case JSON: a number, or an expression string.
"""
param_value(x::Real) = x isa Bool ? throw(ArgumentError("boolean is not a numeric parameter")) :
                       Float64(x)
param_value(x::AbstractString) = parse_param_expr(x)
param_value(x) = throw(ArgumentError("not a numeric parameter: $(repr(x))"))
