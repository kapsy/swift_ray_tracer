

import Foundation

// MARK: Number functions

func assertNaN(_ value: Float) {
    assert(!value.isNaN)
}

func frand48() -> Float {
    return Float(drand48())
}

// MARK: V3

struct V3
{
    var x: Float
    var y: Float
    var z: Float

    // NOTE: (Kapsy) For C union like behavior.
    var r: Float { get { return x } }
    var g: Float { get { return y } }
    var b: Float { get { return z } }

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

// MARK: Camera

struct Camera {

    var origin: V3
    var lowerLeft: V3
    var horiz: V3
    var vert: V3

    var w: V3
    var u: V3
    var v: V3

    var lensRad: Float


    init() {
        self.origin = V3(0)
        self.lowerLeft = V3(0)
        self.horiz = V3(0)
        self.vert = V3(0)

        self.w = V3(0)
        self.u = V3(0)
        self.v = V3(0)

        self.lensRad = Float(0)
    }
}


// MARK: Ray

struct Ray {

    var A: V3
    var B: V3

    // NOTE: (Kapsy) For C union like behavior.
    var orig: V3 { get { return A } }
    var dir: V3 { get { return B } }

    init(_ A: V3, _ B: V3) {
        self.A = A
        self.B = B
    }
}

func getRay(_ c: inout Camera, _ s: Float, _ t: Float) -> Ray {

    // NOTE: (Kapsy) Random in unit disk.
    var rand = V3(0)
    repeat  {
        rand = 2.0*V3(frand48(), frand48(), 0) - V3(1,1,0)
    } while dot(rand, rand) >= 1.0

    let rd = c.lensRad*rand;
    let offset = c.u*rd.x + c.v*rd.y

    let res = Ray(c.origin + offset, c.lowerLeft + s*c.horiz + t*c.vert - c.origin - offset)
    return res
}

func getColor(_ ray: inout Ray, _ depth: Int) -> V3 {
    return V3(1,0,0)
}

func main() {


    let frameRate = Float(25)
    let durationS = Float(5)
    let frameCount = 1//frameRate*durationS

    let testStart = Float(-2)
    let testEnd = Float(2)
    let testDelta = (testEnd - testStart)/Float(frameCount)

    let lookFromTest = testStart

    for f in 0..<frameCount {

        var data = Data()

        let nx = Int(200)
        let ny = Int(100)
        let ns = Int(1)

        // NOTE: (Kapsy) Camera setup stuff.

        let lookFrom = V3(0,1,2)
        let lookAt = V3(0,0.25,-0.2)
        let vup = V3(0,1,0)
        let vfov = Float(50)
        let aspect = Float(nx)/Float(ny)
        let aperture = Float(0)//(0.04)
        let focusDist = length(lookFrom - lookAt)

        let theta = vfov*Float.pi/180
        let halfHeight = tan(theta/2)
        let halfWidth = Float(aspect*halfHeight)

        var cam = Camera()

        cam.origin = lookFrom
        cam.w = unit(lookFrom - lookAt)
        cam.u = unit(cross(vup, cam.w))
        cam.v = cross(cam.w, cam.u)

        cam.lowerLeft = cam.origin - halfWidth*focusDist*cam.u - halfHeight*focusDist*cam.v - focusDist*cam.w
        cam.horiz = 2*halfWidth*focusDist*cam.u
        cam.vert = 2*halfHeight*focusDist*cam.v


        var lensRad: Float

        let header = "P3\n\(nx) \(ny)\n255\n".data(using: .ascii)!
        data.append(header)

        for j in (0..<ny).reversed() {

            for i in 0..<nx {

                var col = V3(0)

                for s in 0..<ns {

                    let u = (Float(i) + frand48())/Float(nx)
                    let v = (Float(j) + frand48())/Float(ny)

                    var r = getRay(&cam, u, v)

                    var depth = 0

                    col += getColor(&r, depth)

                    assertNaN(col.r)
                    assertNaN(col.g)
                    assertNaN(col.b)
                }

                col /= Float(ns)

                let ir = Int(255.99*col.r)
                let ig = Int(255.99*col.g)
                let ib = Int(255.99*col.b)

                assert((ir >= 0) && (ir <= 255))
                assert((ig >= 0) && (ig <= 255))
                assert((ib >= 0) && (ib <= 255))

                let pixel = "\(ir) \(ig) \(ib)\n".data(using: .ascii)!
                data.append(pixel)
            }
        }

        let outPath = "temp/out_\(String(format: "%03d", f + 1)).ppm"

        FileManager.default.createFile(atPath: outPath, contents: data)
    }
}

main()
