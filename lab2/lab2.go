package main

import (
	"fmt"
	"image"
	"image/color"
	"image/gif"
	"log"
	"math"

	"github.com/gonum/matrix/mat64"
	"github.com/llgcode/draw2d/draw2dimg"
)

type Vih struct {
	ind int
	xm0 []float64
	ym0 []float64
	g   []float64
}

func scl(a, b, c, d float64) float64 {
	return a*c + b*d
}

func dist(a, b, c, d float64) float64 {
	return math.Sqrt((a-c)*(a-c) + (b-d)*(b-d))
}

func main() {
	var img gif.GIF

	N := 300
	dt := 0.03
	M := 100
	u0, v0 := 2., 0.
	G0 := 0.
	p := []*Vih{
		{ind: 0},
		//{ind: 49},
		{ind: 99},
	}
	R := 9 * math.Phi / float64(M)

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

	xl0 := make([]float64, M)
	yl0 := make([]float64, M)

	for i := range xl0 {
		xl0[i], yl0[i] = L(i)
	}

	xlc := make([]float64, M-1)
	ylc := make([]float64, M-1)

	for i := range xlc {
		xlc[i] = (xl0[i] + xl0[i+1]) / 2
		ylc[i] = (yl0[i] + yl0[i+1]) / 2
	}

	delta := math.Pow(xl0[1]-xl0[0], 2) + math.Pow(yl0[1]-yl0[0], 2)
	for j := range xlc {
		delta0 := math.Pow(xl0[j+1]-xl0[j], 2) + math.Pow(yl0[j+1]-yl0[j], 2)
		if delta0 < delta {
			delta = delta0
		}
	}
	delta = 0.5 * math.Sqrt(delta)

	norm := func(k int) (float64, float64) {
		d := math.Sqrt(math.Pow(xl0[k+1]-xl0[k], 2) + math.Pow(yl0[k+1]-yl0[k], 2))
		txk := (xl0[k+1] - xl0[k]) / d
		tyk := (yl0[k+1] - yl0[k]) / d

		return -tyk, txk
	}

	Vg := func(x, y, x0, y0 float64) (float64, float64) {
		R0j := math.Max(delta, math.Sqrt(math.Pow(x-x0, 2)+math.Pow(y-y0, 2)))
		return (y0 - y) / (2 * math.Pi * R0j * R0j), (x - x0) / (2 * math.Pi * R0j * R0j)
	}

	Vj := func(x, y float64, j int) (float64, float64) {
		return Vg(x, y, xl0[j], yl0[j])
	}

	// - Time iteration

	for n := 0; n < N; n++ {
		A := mat64.NewDense(M, M, nil)
		b := mat64.NewDense(M, 1, nil)
		G := &mat64.Dense{}

		for k := range xlc {
			nx, ny := norm(k)
			for j := range xl0 {
				u, v := Vj(xlc[k], ylc[k], j)
				A.Set(k, j, scl(u, v, nx, ny))
			}

			rightPart := -scl(u0, v0, nx, ny)

			for _, vih := range p {
				for i, g := range vih.g {
					//u, v := Vj(vih.xm0[i], vih.ym0[i], k)
					u, v := Vg(xlc[k], ylc[k], vih.xm0[i], vih.ym0[i])
					rightPart -= g * scl(u, v, nx, ny)
				}
			}

			b.Set(k, 0, rightPart)
		}

		for j := range xl0 {
			A.Set(M-1, j, 1)
		}
		rightPart := G0
		for _, vih := range p {
			for _, g := range vih.g {
				rightPart -= g
			}
		}
		b.Set(M-1, 0, G0)

		fmt.Println(n)
		//fmt.Println(A)

		if err := G.Solve(A, b); err != nil {
			log.Fatal(err)
		}

		// Answer

		fi := func(x, y float64) float64 {
			var ans float64
			var Gk float64

			for j := range xlc {
				Gk += b.At(j, 0)

				Rj := math.Max(delta, math.Sqrt(math.Pow(x-xlc[j], 2)+math.Pow(y-ylc[j], 2)))

				ans += Gk * ((yl0[j+1]-yl0[j])*(x-xlc[j]) - (xl0[j+1]-xl0[j])*(y-ylc[j])) / (2. * math.Pi * Rj * Rj)
			}

			ans += G0 / (2 * math.Pi) * math.Atan2(y-yl0[M-1], x-xl0[M-1])

			return x*v0 + y*u0 + ans
		}

		V := func(x, y float64) (float64, float64) {
			ansU, ansV := u0, v0

			//fmt.Println("-----------------------------")
			for j := range xl0 {
				uj, vj := Vj(x, y, j)
				ansU += G.At(j, 0) * uj
				ansV += G.At(j, 0) * vj
			}

			for _, vih := range p {
				for i, g := range vih.g {
					u, v := Vg(x, y, vih.xm0[i], vih.ym0[i])
					//u, v := Vg(vih.xm0[i], vih.ym0[i], x, y)
					ansU += g * u
					ansV += g * v
				}
			}

			return ansU, ansV
		}

		// Graphics

		if n%5 == 0 {
			img := image.NewRGBA(image.Rect(0, 0, 1600, 1600))
			for i := 0; i < 1600; i++ {
				for j := 0; j < 1600; j++ {
					img.Set(i, j, color.White)
				}
			}

			gc := draw2dimg.NewGraphicContext(img)

			gc.SetFillColor(color.RGBA{0, 0, 0, 0})
			gc.SetStrokeColor(color.RGBA{0x44, 0x44, 0x44, 0xff})
			gc.SetLineWidth(2)

			for x := 0.; x < 8; x += 0.15 {
				for y := 0.; y < 8; y += 0.15 {
					fi(x, y)
					u, v := V(x, y)
					//d := math.Sqrt(u*u + v*v)

					gc.MoveTo(x*200, 1600-y*200)
					gc.LineTo(x*200+u*5, 1600-(y*200+v*5))

					//fmt.Printf("%5.2f %5.2f | %6.3f %6.3f - %5.3f\n", x, y, u, v, d)
				}
			}

			gc.FillStroke()

			gc.SetFillColor(color.RGBA{0, 0, 0, 0})
			gc.SetStrokeColor(color.RGBA{0xff, 0, 0, 0xff})
			gc.SetLineWidth(2)

			for i := 0; i < M; i++ {
				x, y := L(i)

				gc.MoveTo(x*200-3, 1600-(y*200-3))
				gc.LineTo(x*200+3, 1600-(y*200+3))
				gc.MoveTo(x*200-3, 1600-(y*200+3))
				gc.LineTo(x*200+3, 1600-(y*200-3))
			}

			gc.FillStroke()

			gc.SetStrokeColor(color.RGBA{0, 0, 0xff, 0xff})
			gc.SetLineWidth(2)

			for _, vih := range p {
				for i := range vih.g {
					x, y := vih.xm0[i], vih.ym0[i]

					gc.MoveTo(x*200-3, 1600-(y*200-3))
					gc.LineTo(x*200+3, 1600-(y*200+3))
					gc.MoveTo(x*200-3, 1600-(y*200+3))
					gc.LineTo(x*200+3, 1600-(y*200-3))
				}
			}

			gc.Close()
			gc.FillStroke()

			draw2dimg.SaveToPngFile(fmt.Sprintf("img%d.png", n), img)
		}
		// - End of Graphics

		// -- Adding end points
		for _, vih := range p {
			vih.xm0 = append(vih.xm0, xl0[vih.ind])
			vih.ym0 = append(vih.ym0, yl0[vih.ind])
			vih.g = append(vih.g, G.At(vih.ind, 0))
		}

		// -- Calculating new positions
		for _, vih := range p {
			for i := range vih.g {
				u, v := V(vih.xm0[i], vih.ym0[i])
				vih.xm0[i] = vih.xm0[i] + u*dt
				vih.ym0[i] = vih.ym0[i] + v*dt
			}
		}

		// -- Moving away from line

		for _, vih := range p {
			for i := range vih.g {
				for j := range xlc {
					x0, y0 := xlc[j], ylc[j]
					x, y := vih.xm0[i], vih.ym0[i]

					d := dist(x0, y0, x, y)
					if d < R {
						dx := x - x0
						dy := y - y0

						n := math.Sqrt(dx*dx + dy*dy)
						x = x0 + (x-x0)*R/n
						y = y0 + (y-y0)*R/n
					}

					vih.xm0[i], vih.ym0[i] = x, y
				}
			}
		}
	}
}
