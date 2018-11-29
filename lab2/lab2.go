package main

import (
	"flag"
	"image"
	"image/color"
	"image/gif"
	"log"
	"math"
	"os"
	"sync"
	"time"

	"github.com/gonum/matrix/mat64"
)

var (
	threads    = flag.Int("threads", 12, "Number of threads for computing.")
	iterations = flag.Int("iterations", 1500, "Number of iterations.")
	dt         = flag.Float64("dt", 0.1, "Time range for one iteration.")
)

type Vih struct {
	ind    int
	xm0    []float64
	ym0    []float64
	preXM0 []float64
	preYM0 []float64
	u      []float64
	v      []float64
	g      []float64
}

func scl(a, b, c, d float64) float64 {
	return a*c + b*d
}

func dist(a, b, c, d float64) float64 {
	return math.Sqrt((a-c)*(a-c) + (b-d)*(b-d))
}

func main() {
	flag.Parse()

	var img gif.GIF

	var (
		matrixTime   time.Duration
		graphicsTime time.Duration
		pointsTime   time.Duration
		movingTime   time.Duration
	)

	N := *iterations
	dt := *dt
	M := 51
	u0, v0 := 1., 0.
	G0 := 0.
	p := []*Vih{
		{ind: 0},
		{ind: M - 1},
	}
	R := math.Phi / float64(M)

	L := func(i int) (float64, float64) {
		mm := (M + 1) / 2
		var (
			angle float64
			dy    float64
		)
		if i < mm {
			angle = float64(i)/float64(mm-1)*math.Pi*1.1 + math.Pi*0.4
			dy = 1
		} else {
			angle = float64(i-mm+1)/float64(M-mm)*math.Pi*1.1 + math.Pi*0.5
			dy = -1
		}
		return math.Cos(angle), math.Sin(angle) + dy
	}

	xl0 := make([]float64, M)
	yl0 := make([]float64, M)

	for i := range xl0 {
		xl0[i], yl0[i] = L(i)
		log.Printf("%6.2f %6.2f\n", xl0[i], yl0[i])
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
		dx := xl0[k+1] - xl0[k]
		dy := yl0[k+1] - yl0[k]
		d := math.Sqrt(dx*dx + dy*dy)

		return -dy / d, dx / d
	}

	//Vg := func(x, y, x0, y0 float64) (float64, float64) {
	//	R0j2 := math.Max(delta*delta, math.Pow(x-x0, 2)+math.Pow(y-y0, 2))
	//	return (y0 - y) / (2 * math.Pi * R0j2), (x - x0) / (2 * math.Pi * R0j2)
	//}

	delta2 := delta * delta
	Vg := func(x, y, x0, y0 float64) (float64, float64) {
		R0j2 := math.Max(delta2, (x-x0)*(x-x0)+(y-y0)*(y-y0))
		R0j2 = 1 / (2 * math.Pi * R0j2)
		return (y0 - y) * R0j2, (x - x0) * R0j2
	}

	Vj := func(x, y float64, j int) (float64, float64) {
		return Vg(x, y, xl0[j], yl0[j])
	}

	// - Time iteration

	A := mat64.NewDense(M, M, nil)
	b := mat64.NewDense(M, 1, nil)
	G := &mat64.Dense{}
	InvA := &mat64.Dense{}

	for k := range xlc {
		nx, ny := norm(k)
		for j := range xl0 {
			u, v := Vj(xlc[k], ylc[k], j)
			A.Set(k, j, scl(u, v, nx, ny))
		}
	}
	for j := range xl0 {
		A.Set(M-1, j, 1)
	}

	if err := InvA.Inverse(A); err != nil {
		log.Fatal(err)
	}

	type RemoteVTask struct {
		i, j int
		V    func(float64, float64) (float64, float64)
	}

	wg := sync.WaitGroup{}
	var tasks = make(chan int, 100)
	var remoteVTasks = make(chan RemoteVTask, 100)
	for i := 0; i < *threads; i++ {
		go func() {
			for k := range tasks {
				func() {
					defer wg.Done()
					nx, ny := norm(k)
					rightPart := -scl(u0, v0, nx, ny)

					for _, vih := range p {
						for i, g := range vih.g {
							u, v := Vg(xlc[k], ylc[k], vih.xm0[i], vih.ym0[i])
							rightPart -= g * scl(u, v, nx, ny)
						}
					}

					b.Set(k, 0, rightPart)
				}()
			}
		}()

		go func() {
			for t := range remoteVTasks {
				vih := p[t.j]
				vih.u[t.i], vih.v[t.i] = t.V(vih.xm0[t.i], vih.ym0[t.i])
				wg.Done()
			}
		}()
	}

	tt0 := time.Now()
	t00 := time.Now()
	totalFrameTime := 0.
	lastFrameTime := 0.
	for n := 0; n < N; n++ {
		log.Println(n, time.Now().Sub(t00))
		t00 = time.Now()

		// Matrix time
		t0 := time.Now()
		for k := range xlc {
			wg.Add(1)
			tasks <- k
		}
		rightPart := G0
		for _, vih := range p {
			for _, g := range vih.g {
				rightPart -= g
			}
		}
		b.Set(M-1, 0, rightPart)
		wg.Wait()

		G.Mul(InvA, b)
		matrixTime += time.Now().Sub(t0)
		// End of matrix time

		g := G.RawMatrix().Data

		V := func(x, y float64) (float64, float64) {
			ansU, ansV := u0, v0

			for j := range xl0 {
				uj, vj := Vj(x, y, j)
				ansU += g[j] * uj
				ansV += g[j] * vj
			}

			for _, vih := range p {
				for i, g := range vih.g {
					u, v := Vg(x, y, vih.xm0[i], vih.ym0[i])
					ansU += g * u
					ansV += g * v
				}
			}

			return ansU, ansV
		}

		// Graphics
		if dt := int((totalFrameTime - lastFrameTime) * 100); dt > 5 {
			t0 = time.Now()
			pal := PrintVec(V, L, p, M, u0, v0)
			img.Image = append(img.Image, pal)

			img.Delay = append(img.Delay, dt)
			lastFrameTime += float64(dt) / 100.

			graphicsTime += time.Now().Sub(t0)
		}
		// - End of Graphics

		t0 = time.Now()
		// -- Adding end points
		for _, vih := range p {
			vih.xm0 = append(vih.xm0, xl0[vih.ind])
			vih.ym0 = append(vih.ym0, yl0[vih.ind])
			vih.preXM0 = append(vih.preXM0, xl0[vih.ind])
			vih.preYM0 = append(vih.preYM0, yl0[vih.ind])
			vih.u = append(vih.u, 0)
			vih.v = append(vih.v, 0)
			vih.g = append(vih.g, G.At(vih.ind, 0))
		}

		// -- Calculating new positions
		for j, vih := range p {
			for i := range vih.g {
				wg.Add(1)
				remoteVTasks <- RemoteVTask{
					i: i,
					j: j,
					V: V,
				}
			}
		}
		wg.Wait()
		var maxW float64
		for _, vih := range p {
			for i := range vih.g {
				w := vih.u[i]*vih.u[i] + vih.v[i]*vih.v[i]
				if w > maxW {
					maxW = w
				}
			}
		}
		mdt := delta / maxW
		totalFrameTime += mdt
		if mdt > dt {
			mdt = dt
		}
		for _, vih := range p {
			for i := range vih.g {
				vih.preXM0[i] = vih.xm0[i]
				vih.preYM0[i] = vih.ym0[i]
				vih.xm0[i] += vih.u[i] * mdt
				vih.ym0[i] += vih.v[i] * mdt
			}
		}

		pointsTime += time.Now().Sub(t0)

		// -- Moving away from line
		t0 = time.Now()

		// Old
		//for _, vih := range p {
		//	for i := range vih.g {
		//		minD := R * 2
		//		minJ := -1
		//		for j := range xlc {
		//			x0, y0 := xlc[j], ylc[j]
		//			x, y := vih.xm0[i], vih.ym0[i]
		//
		//			d := dist(x0, y0, x, y)
		//			if d < minD {
		//				minD = d
		//				minJ = j
		//			}
		//		}
		//		if minD < R {
		//			x0, y0 := xlc[minJ], ylc[minJ]
		//			x, y := vih.xm0[i], vih.ym0[i]
		//
		//			dx := x - x0
		//			dy := y - y0
		//
		//			n := math.Sqrt(dx*dx + dy*dy)
		//			x = x0 + dx*R/n
		//			y = y0 + dy*R/n
		//			vih.xm0[i], vih.ym0[i] = x, y
		//		}
		//	}
		//}

		// New
		t0 = time.Now()
		for _, vih := range p {
			for i := range vih.g {
				minD := R * 2
				minJ := -1
				for j := range xlc {
					d := dist(xlc[j], ylc[j], vih.xm0[i], vih.ym0[i])
					if d < minD {
						minD = d
						minJ = j
					}
				}
				if minD < R {
					xc, yc := xlc[minJ], ylc[minJ]
					x0, y0 := vih.preXM0[i], vih.preYM0[i]
					x1, y1 := vih.xm0[i], vih.ym0[i]

					a, b, c := x1-x0, y1-y0, (x0-x1)*xc+(y0-y1)*yc
					k := 1. / math.Sqrt(a*a+b*b)

					d0 := k * (a*x0 + b*y0 + c)
					d1 := k * (a*x1 + b*y1 + c)

					if d0*d1 < 0 {
						x0, y0 = x1, x1
					}
					dx := x0 - xc
					dy := y0 - yc

					n := math.Sqrt(dx*dx + dy*dy)
					vih.xm0[i], vih.ym0[i] = xc+dx*R/n, yc+dy*R/n
				}
			}
		}

		movingTime += time.Now().Sub(t0)
	}
	totalTime := time.Now().Sub(tt0)

	log.Printf("matrixTime: %5.2f%% %v\n", float64(matrixTime*100)/float64(totalTime), matrixTime)
	log.Printf("graphicsTime: %5.2f%% %v\n", float64(graphicsTime*100)/float64(totalTime), graphicsTime)
	log.Printf("pointsTime: %5.2f%% %v\n", float64(pointsTime*100)/float64(totalTime), pointsTime)
	log.Printf("movingTime: %5.2f%% %v\n", float64(movingTime*100)/float64(totalTime), movingTime)
	log.Println("--------------------------")
	log.Printf("totalTime: %v\n", totalTime)

	f, _ := os.Create("img.gif")
	if err := gif.EncodeAll(f, &img); err != nil {
		log.Print(err)
	}
}

const (
	sizeX = 1200
	sizeY = 400
	k     = 200 / 2
	//minX, maxX = -2.5, 2.5
	//minY, maxY = -2.5, 2.5
	minX, maxX = -2., 22.
	minY, maxY = -4., 4.
)

var (
	colorWhite = color.White
	colorRed   = color.RGBA{255, 0, 0, 255}
	colorGrey  = color.RGBA{0x44, 0x44, 0x44, 0xff}
	colorBlue  = color.RGBA{0, 0, 0xff, 0xff}
)

func convCoord(x0, y0 float64) (x, y int) {
	x = int((x0 - minX) * float64(sizeX) / (maxX - minX))
	y = sizeY - int((y0-minY)*float64(sizeY)/(maxY-minY))
	return
}

func PrintVec(V func(x, y float64) (float64, float64), L func(int) (float64, float64), p []*Vih, M int, u0, v0 float64) *image.Paletted {
	img := image.NewPaletted(image.Rect(0, 0, sizeX, sizeY), color.Palette{colorWhite, colorGrey, colorRed, colorBlue})

	//d0 := math.Sqrt(u0*u0 + v0*v0)
	//for x := minX; x < maxX; x += 0.25 {
	//	for y := minY; y < maxY; y += 0.25 {
	//		u, v := V(x, y)
	//		d := math.Sqrt(u*u + v*v)
	//		if d > 3*d0 {
	//			u *= 3 / d
	//			v *= 3 / d
	//		} else {
	//			u /= d0
	//			v /= d0
	//		}
	//		u *= 10
	//		v *= 10
	//
	//		x0, y0 := convCoord(x, y)
	//		x1, y1 := convCoord(x+u/k, y+v/k)
	//		Line(img, x0, y0, x1, y1, 1)
	//	}
	//}

	for i := 0; i < M; i++ {
		x, y := L(i)

		x0, y0 := convCoord(x, y)
		Tic(img, x0, y0, 2)
	}

	for _, vih := range p {
		for i := range vih.g {
			x, y := vih.xm0[i], vih.ym0[i]

			x0, y0 := convCoord(x, y)
			Tic(img, x0, y0, 3)
		}
	}

	return img
}

func Line(img *image.Paletted, x0, y0, x, y int, ind uint8) {
	d1 := x - x0
	d2 := y - y0
	if d1 < 0 {
		d1 = -d1
	}
	if d2 < 0 {
		d2 = -d2
	}

	d := d1
	if d < d2 {
		d = d2
	}

	for i := 0; i < d; i++ {
		x1 := x0 + i*(x-x0)/d
		y1 := y0 + i*(y-y0)/d
		img.SetColorIndex(x1, y1, ind)
	}
}

func Tic(img *image.Paletted, x, y int, ind uint8) {
	Line(img, x-3, y-3, x+3, y+3, ind)
	Line(img, x-3, y+3, x+3, y-3, ind)
}
