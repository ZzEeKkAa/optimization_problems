package main

import (
	"fmt"
	"image"
	"image/color"
	"log"
	"math"

	"github.com/gonum/matrix/mat64"
	"github.com/llgcode/draw2d/draw2dimg"
)

func scl(a, b, c, d float64) float64 {
	return a*c + b*d
}

var u0, v0 float64

func main() {
	M := 100
	//u0, v0 = 0.2, 0.1
	u0, v0 = -0.2, 0.0
	G0 := 0.

	L := func(i int) (float64, float64) {
		mm := (M + 1) / 2
		var (
			angle float64
			dy    float64
		)
		if i < mm {
			angle = float64(i)/float64(mm-1)*math.Pi*1.3 + math.Pi*0.2
			dy = 4
		} else {
			angle = float64(i-mm+1)/float64(M-mm)*math.Pi*1.3 + math.Pi*0.5
			dy = 2
		}
		return math.Cos(angle) + 4, math.Sin(angle) + dy + 1
	}

	x0 := make([]float64, M)
	y0 := make([]float64, M)

	for i := range x0 {
		x0[i], y0[i] = L(i)
	}

	xc := make([]float64, M-1)
	yc := make([]float64, M-1)

	for i := range xc {
		xc[i] = (x0[i] + x0[i+1]) / 2
		yc[i] = (y0[i] + y0[i+1]) / 2
	}

	delta := math.Pow(x0[1]-x0[0], 2) + math.Pow(y0[1]-y0[0], 2)
	for j := range xc {
		delta0 := math.Pow(x0[j+1]-x0[j], 2) + math.Pow(y0[j+1]-y0[j], 2)
		if delta0 < delta {
			delta = delta0
		}
	}
	delta = 0.5 * math.Sqrt(delta)

	n := func(k int) (float64, float64) {
		d := math.Sqrt(math.Pow(x0[k+1]-x0[k], 2) + math.Pow(y0[k+1]-y0[k], 2))
		txk := (x0[k+1] - x0[k]) / d
		tyk := (y0[k+1] - y0[k]) / d

		return -tyk, txk
	}

	Vj := func(x, y float64, j int) (float64, float64) {
		R0j := math.Max(delta, math.Sqrt(math.Pow(x-x0[j], 2)+math.Pow(y-y0[j], 2)))
		return (y0[j] - y) / (2 * math.Pi * R0j * R0j), (x - x0[j]) / (2 * math.Pi * R0j * R0j)
	}

	A := mat64.NewDense(M, M, nil)
	b := mat64.NewDense(M, 1, nil)
	G := &mat64.Dense{}

	for k := range xc {
		nx, ny := n(k)
		for j := range x0 {
			u, v := Vj(xc[k], yc[k], j)
			A.Set(k, j, scl(u, v, nx, ny))
		}
		b.Set(k, 0, -scl(u0, v0, nx, ny))
	}

	for j := range x0 {
		A.Set(M-1, j, 1)
	}
	b.Set(M-1, 0, G0)

	if err := G.Solve(A, b); err != nil {
		log.Fatal(err)
	}

	// Answer

	fi := func(x, y float64) float64 {
		var ans float64
		var Gk float64

		for j := range xc {
			Gk += b.At(j, 0)

			Rj := math.Max(delta, math.Sqrt(math.Pow(x-xc[j], 2)+math.Pow(y-yc[j], 2)))

			ans += Gk * ((y0[j+1]-y0[j])*(x-xc[j]) - (x0[j+1]-x0[j])*(y-yc[j])) / (2. * math.Pi * Rj * Rj)
		}

		ans += G0 / (2 * math.Pi) * math.Atan2(y-y0[M-1], x-x0[M-1])

		return x*v0 + y*u0 + ans
	}

	psi := func(x, y float64) float64 {
		var ans float64
		var Gk float64

		for j := range xc {
			Gk = b.At(j, 0)

			Rj := math.Max(delta, math.Sqrt(math.Pow(x-xc[j], 2)+math.Pow(y-yc[j], 2)))

			ans += Gk / (2. * math.Pi) * math.Log(Rj)
		}

		return -x*v0 + y*u0 - ans
	}
	F := func(x, y float64) (float64, float64) {
		return fi(x, y) / 8, psi(x, y) / 8
	}

	V := func(x, y float64) (float64, float64) {
		ansU, ansV := u0, v0

		for j := range x0 {
			uj, vj := Vj(x, y, j)
			ansU += G.At(j, 0) * uj
			ansV += G.At(j, 0) * vj
		}

		return ansU, ansV
	}
	fi(0, 0)
	V(0, 0)
	F(0, 0)

	// Graphics
	PrintVec(V, L, M)
	//PrintVecLines(V, L, M, 1.6, 0.8)
	//PrintVecLines(F, L, M, 0.8, 3)
	//PrintPressure(V, L, M)
}

const (
	size = 800
	k    = 200 / 2
)

func PrintVec(V func(x, y float64) (float64, float64), L func(int) (float64, float64), M int) {
	img := image.NewRGBA(image.Rect(0, 0, size, size))
	for i := 0; i < size; i++ {
		for j := 0; j < size; j++ {
			img.Set(i, j, color.White)
		}
	}

	gc := draw2dimg.NewGraphicContext(img)

	gc.SetFillColor(color.RGBA{0, 0, 0, 0})
	gc.SetStrokeColor(color.RGBA{0x44, 0x44, 0x44, 0xff})
	gc.SetLineWidth(2)

	d0 := math.Sqrt(u0*u0 + v0*v0)
	for x := 0.; x < 8; x += 0.25 {
		for y := 0.; y < 8; y += 0.25 {
			u, v := V(x, y)
			d := math.Sqrt(u*u + v*v)
			if d > 3*d0 {
				u *= 3 / d
				v *= 3 / d
			} else {
				u /= d0
				v /= d0
			}
			u *= 10
			v *= 10
			gc.MoveTo(x*k, size-y*k)
			gc.LineTo(x*k+u, size-(y*k+v))
		}
	}

	gc.FillStroke()

	gc.SetFillColor(color.RGBA{0, 0, 0, 0})
	gc.SetStrokeColor(color.RGBA{0xff, 0, 0, 0xff})
	gc.SetLineWidth(2)

	for i := 0; i < M; i++ {
		x, y := L(i)

		gc.MoveTo(x*k-3, size-(y*k-3))
		gc.LineTo(x*k+3, size-(y*k+3))
		gc.MoveTo(x*k-3, size-(y*k+3))
		gc.LineTo(x*k+3, size-(y*k-3))
	}

	gc.Close()
	gc.FillStroke()

	draw2dimg.SaveToPngFile("img.png", img)
}

func PrintPressure(V func(x, y float64) (float64, float64), L func(int) (float64, float64), M int) {
	d0 := math.Sqrt(u0*u0 + v0*v0)
	img := image.NewRGBA(image.Rect(0, 0, size, size))
	for i := 0; i < size; i++ {
		for j := 0; j < size; j++ {
			x, y := float64(i)/k, float64(j)/k
			u, v := V(x, y)
			d := math.Sqrt(u*u + v*v)
			if d > 3*d0 {
				d = 1
			} else {
				d /= 3 * d0
			}
			img.Set(i, j, color.RGBA{R: uint8(255 * d), B: uint8(255 * (1 - d)), A: 255})
		}
	}

	gc := draw2dimg.NewGraphicContext(img)

	gc.SetFillColor(color.RGBA{0, 0, 0, 0})
	gc.SetStrokeColor(color.RGBA{0xff, 0, 0, 0xff})
	gc.SetLineWidth(2)

	for i := 0; i < M; i++ {
		x, y := L(i)

		gc.MoveTo(x*k-3, size-(y*k-3))
		gc.LineTo(x*k+3, size-(y*k+3))
		gc.MoveTo(x*k-3, size-(y*k+3))
		gc.LineTo(x*k+3, size-(y*k-3))
	}

	gc.Close()
	gc.FillStroke()

	draw2dimg.SaveToPngFile("img.png", img)
}

func PrintVecLines(V func(x, y float64) (float64, float64), L func(int) (float64, float64), M int, dx, dy float64) {
	img := image.NewRGBA(image.Rect(0, 0, size, size))
	for i := 0; i < size; i++ {
		for j := 0; j < size; j++ {
			img.Set(i, j, color.White)
		}
	}

	gc := draw2dimg.NewGraphicContext(img)

	gc.SetFillColor(color.RGBA{0, 0, 0, 0})
	gc.SetStrokeColor(color.RGBA{0x44, 0x44, 0x44, 0xff})
	gc.SetLineWidth(2)

	buildLine := func(x, y float64) {
		for i := 0; i < 200; i++ {
			gc.MoveTo(x*k, size-y*k)
			u, v := V(x, y)
			fmt.Println(u, v)
			x, y = x+u, y+v
			gc.LineTo(x*k, size-(y*k))
			if x > 8 || y > 8 || x < 0 || y < 0 {
				return
			}
		}
	}

	for x := 0.; x < 8; x = x + dx {
		buildLine(x, 0.)
	}
	for y := 0.; y < 8; y = y + dy {
		buildLine(0., y)
	}

	gc.FillStroke()

	gc.SetFillColor(color.RGBA{0, 0, 0, 0})
	gc.SetStrokeColor(color.RGBA{0xff, 0, 0, 0xff})
	gc.SetLineWidth(2)

	for i := 0; i < M; i++ {
		x, y := L(i)

		gc.MoveTo(x*k-3, size-(y*k-3))
		gc.LineTo(x*k+3, size-(y*k+3))
		gc.MoveTo(x*k-3, size-(y*k+3))
		gc.LineTo(x*k+3, size-(y*k-3))
	}

	gc.Close()
	gc.FillStroke()

	draw2dimg.SaveToPngFile("img.png", img)
}
