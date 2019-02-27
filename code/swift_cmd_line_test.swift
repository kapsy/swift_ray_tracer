import Foundation

func assertNaN(_ value: Float) {
    assert(!value.isNaN)
}

func frand48() -> Float {
    return Float(drand48())
}

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


func /(a: V3, b: Float) -> V3 {
    let res = V3(a.x/b, a.y/b, a.z/b)
    return res
}

func /=(a: inout V3, b: Float) {
    a = a/b
}

func +(a: V3, b: V3) -> V3 {
    let res = V3(a.x + b.x, a.y + b.y, a.z + b.z)
    return res
}

func +=(a: inout V3, b: V3) {
    a = a + b
}

func main() {

    let frameCount = 1

    for f in 0..<frameCount {

        var data = Data()

        let nx = Int(200)
        let ny = Int(100)
        let ns = Int(1)

        let header = "P3\n\(nx) \(ny)\n255\n".data(using: .ascii)!
        data.append(header)

        for j in (0..<ny).reversed()
        {
            for i in 0..<nx
            {

                var col = V3(0)

                for s in 0..<ns {

                    col += V3(Float(i)/Float(nx),
                              Float(j)/Float(ny),
                              Float(0.2))

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
