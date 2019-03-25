import Foundation

// MARK: Number functions

func assertNaN(_ value: Float) {
    assert(!value.isNaN)
}

func filterNaN(_ value: Float) -> Float {
    return value.isNaN ? 0.0 : value
}

func drand48f() -> Float {
    return Float(drand48())
}

func clamp01(_ a: Float) -> Float {
    if a > 1.0 {
        return 1.0
    } else if a < 0.0 {
        return 0.0
    }
    return a
}

// MARK: V3

struct V3 {
    var x = Float(0)
    var y = Float(0)
    var z = Float(0)

    // NOTE: (Kapsy) For C union like behavior.
    var r: Float { get { return x } set(a) { self.x = a } }
    var g: Float { get { return y } set(a) { self.y = a } }
    var b: Float { get { return z } set(a) { self.z = a } }


    init() {}

    init(_ val: Float) {
        self.x = val
        self.y = val
        self.z = val
    }

    init(_ x: Float, _ y: Float, _ z: Float) {
        self.x = x
        self.y = y
        self.z = z
    }
}

func +(a: V3, b: V3) -> V3 {
    let res = V3(a.x + b.x, a.y + b.y, a.z + b.z)
    return res
}

func +=(a: inout V3, b: V3) {
    a = a + b
}

func -(a: V3, b: V3) -> V3 {
    let res = V3(a.x - b.x, a.y - b.y, a.z - b.z)
    return res
}

func -=(a: inout V3, b: V3) {
    a = a - b
}

func *(a: V3, b: V3) -> V3 {
    let res = V3(a.x*b.x, a.y*b.y, a.z*b.z)
    return res
}

func *=(a: inout V3, b: V3) {
    a = a * b
}

func /(a: V3, b: Float) -> V3 {
    let res = V3(a.x/b, a.y/b, a.z/b)
    return res
}

func /=(a: inout V3, b: Float) {
    a = a/b
}

func *(a: Float, b: V3) -> V3 {
    let res = V3(a*b.x, a*b.y, a*b.z)
    return res
}

func *(a: V3, b: Float) -> V3 {
    let res = b*a
    return res
}

prefix func -(a: V3) -> V3 {
    let res = V3(-a.x, -a.y, -a.z)
    return res
}

func dot(_ a: V3, _ b: V3) -> Float {
    let res = a.x*b.x + a.y*b.y + a.z*b.z
    return res
}

func cross(_ a: V3, _ b: V3) -> V3 {
    let res = V3(a.y*b.z - a.z*b.y,
                 a.z*b.x - a.x*b.z,
                 a.x*b.y - a.y*b.x);
    return res
}

func length(_ a: V3) -> Float {
    let res = sqrt(dot(a, a))
    return res
}

func unit(_ a: V3) -> V3 {
    let res = a/length(a)
    return res
}

func squaredLen(_ v: V3) -> Float {
    return dot(v, v)
}

func randomInUnitSphere() -> V3 {
    var v = V3(0)

    repeat {
        v = 2.0*V3(drand48f(), drand48f(), drand48f()) - V3(1)
    } while squaredLen(v) >= 1.0

    return (v);
}
