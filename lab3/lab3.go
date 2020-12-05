package main

import (
	"flag"
	"fmt"
	"image"
	"image/color"
	"image/gif"
	"log"
	"math"
	"os"
	"sync"
	"time"

	"golang.org/x/image/font"
	"golang.org/x/image/font/basicfont"
	"golang.org/x/image/math/fixed"

	"github.com/gonum/matrix/mat64"
)

var (
	threads         = flag.Int("threads", 12, "Number of threads for computing.")
	klyaksaPoints   = flag.Int("klyaksa_points", 20, "Number of threads for computing.")
	iterations      = flag.Int("iterations", 2200, "Number of iterations.")
	iterationsStart = flag.Int("iterations_start", 1400, "Iterations from what animation is writing to gif.")
	dt              = flag.Float64("dt", 1, "Time range for one iteration.")
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

type Klyaksa struct {
	x              []float64
	y              []float64
	u              []float64
	v              []float64
	ax, bx, cx, dx []float64
	ay, by, cy, dy []float64
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
		klyaksaTime  time.Duration
		splinesTime  time.Duration
		movingTime   time.Duration
	)

	kl := Klyaksa{}
	for phi := 0.; phi < 2*math.Pi; phi += 2 * math.Pi / float64(*klyaksaPoints) {
		fmt.Println(phi)
		kl.x = append(kl.x, -5+2*math.Cos(phi))
		kl.y = append(kl.y, 2*math.Sin(phi))
		kl.u = append(kl.u, 0)
		kl.v = append(kl.v, 0)
	}

	N := *iterations
	NStart := *iterationsStart
	dt := *dt
	M := 21
	u0, v0 := 1., 0.
	G0 := 0.
	p := []*Vih{
		{ind: 0},
		{ind: M - 1},
	}
	//R := math.Phi / float64(M)

	alpha := math.Pi * 0.2
	betta := math.Pi * 1.1
	dY := math.Cos(alpha)

	L := func(i int) (float64, float64) {
		mm := (M + 1) / 2
		var (
			angle float64
			dy    float64
		)
		if i < mm {
			angle = float64(i)/float64(mm-1)*(betta-alpha) + (1.5*math.Pi - betta)
			dy = dY
		} else {
			angle = float64(i-mm+1)/float64(M-mm)*(betta-alpha) + (math.Pi*0.5 + alpha)
			dy = -dY
		}
		return math.Cos(angle), math.Sin(angle) + dy
	}

	drawSolidLetter := func(pal *image.Paletted) {
		N := 200
		nn := (N + 1) / 2
		for i := 0; i < N; i++ {
			var (
				angle float64
				dy    float64
			)
			if i < nn {
				angle = float64(i)/float64(nn-1)*(betta-alpha) + (1.5*math.Pi - betta)
				dy = dY
			} else {
				angle = float64(i-nn+1)/float64(N-nn)*(betta-alpha) + (math.Pi*0.5 + alpha)
				dy = -dY
			}
			x, y := convCoord(math.Cos(angle), math.Sin(angle)+dy)
			pal.SetColorIndex(x, y, 1)
			pal.SetColorIndex(x+1, y, 1)
			pal.SetColorIndex(x-1, y, 1)
			pal.SetColorIndex(x, y+1, 1)
			pal.SetColorIndex(x, y-1, 1)
			pal.SetColorIndex(x+1, y+1, 1)
			pal.SetColorIndex(x-1, y-1, 1)
			pal.SetColorIndex(x-1, y+1, 1)
			pal.SetColorIndex(x+1, y-1, 1)
		}
	}

	//moveAway := func(x1, y1, x2, y2 float64) (float64, float64) {
	//	a := x1 - x2
	//	b := y1 - y2
	//	c := a * a
	//	d := b * b
	//	e := c + d
	//	f := y2*x1 - x2*y1
	//	g := f * f
	//	h := math.Sqrt(c * (-g + c + d))
	//
	//	x3 := (-f*b + h) / e
	//	y3 := (h*b + c*f) / (e * a)
	//	x4 := (-f*b - h) / e
	//	y4 := (-h*b + c*f) / (e * a)
	//
	//	fmt.Println("| ", x3, y3)
	//	fmt.Println("| ", x4, y4)
	//
	//	return x3, y3
	//}

	//x1, y1, x2, y2 := math.Cos(math.Pi/4)-0.0001+100, math.Sin(math.Pi/4)-0.0001+100, math.Cos(math.Pi/4)+0.0001+100, math.Sin(math.Pi/4)+0.0001+100
	//fmt.Println(x1, y1, x2, y2)
	//moveAway(x1, y1, x2, y2)
	//os.Exit(0)

	moveAway := func(x1, y1, x2, y2 float64) (float64, float64) {
		y1 -= dY
		y2 -= dY
		d1 := x1*x1 + y1*y1
		d2 := x2*x2 + y2*y2
		if d1 < 1 && d2 > 1 || d1 > 1 && d2 < 1 {
			a := (math.Atan2(y1, x1) + math.Atan2(y2, x2)) / 2
			if a < 0 {
				a += math.Pi * 2
			}
			a = 1.5*math.Pi - a
			if alpha < a && a < betta {
				k := d1 / d2
				x2 *= k
				y2 *= k
				y2 += dY
				return x2, y2
			}
		}
		y1 += 2 * dY
		y2 += 2 * dY
		d1 = x1*x1 + y1*y1
		d2 = x2*x2 + y2*y2
		if d1 < 1 && d2 > 1 || d1 > 1 && d2 < 1 {
			a := (math.Atan2(y1, x1) + math.Atan2(y2, x2)) / 2
			if a < 0 {
				a += math.Pi * 2
			}
			a -= 0.5 * math.Pi
			if alpha < a && a < betta {
				k := d1 / d2
				x2 *= k
				y2 *= k
				y2 += dY
				return x2, y2
			}
		}
		y2 -= dY
		return x2, y2
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

	type RemoteKlTask struct {
		i int
		V func(float64, float64) (float64, float64)
	}

	wg := sync.WaitGroup{}
	var tasks = make(chan int, 100)
	var remoteVTasks = make(chan RemoteVTask, 100)
	var remoteKlTasks = make(chan RemoteKlTask, 100)
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

		go func() {
			for t := range remoteKlTasks {
				kl.u[t.i], kl.v[t.i] = t.V(kl.x[t.i], kl.y[t.i])
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

		t0 = time.Now()
		kl.BuildSplines()
		splinesTime += time.Now().Sub(t0)

		// Graphics
		if n >= NStart {
			if n == NStart {
				totalFrameTime = lastFrameTime
			}
			if dt := int((totalFrameTime - lastFrameTime) * 100); dt > 5 {
				t0 = time.Now()
				pal := PrintVec(V, L, p, M, u0, v0)
				addLabel(pal, sizeX-280, 20, "Yevhenii Havrylko, CM-2, September 2019")
				addLabel(pal, sizeX-320, sizeY-20, fmt.Sprintf("dt = %5.2f, t = %7.2fs, iter = %7d", float64(dt)/100, totalFrameTime, n))
				kl.AddSplineToImage(pal)
				img.Image = append(img.Image, pal)

				img.Delay = append(img.Delay, dt)
				lastFrameTime += float64(dt) / 100.

				drawSolidLetter(pal)

				graphicsTime += time.Now().Sub(t0)
			}
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

		pointsTime += time.Now().Sub(t0)
		t0 = time.Now()

		for i := range kl.x {
			wg.Add(1)
			remoteKlTasks <- RemoteKlTask{
				i: i,
				V: V,
			}
		}

		wg.Wait()
		klyaksaTime += time.Now().Sub(t0)
		t0 = time.Now()

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
		t0 = time.Now()

		if n > 1500 {
			for i := range kl.x {
				kl.x[i] += kl.u[i] * mdt
				kl.y[i] += kl.v[i] * mdt
			}
			kl.AddPoints(2 * 2 * math.Pi / float64(*klyaksaPoints))
		}
		klyaksaTime += time.Now().Sub(t0)

		// -- Moving away from line
		t0 = time.Now()

		// New
		t0 = time.Now()
		for _, vih := range p {
			for i := range vih.g {
				vih.xm0[i], vih.ym0[i] = moveAway(vih.preXM0[i], vih.preYM0[i], vih.xm0[i], vih.ym0[i])
				//minD := R * 2
				//minJ := -1
				//for j := range xlc {
				//	d := dist(xlc[j], ylc[j], vih.xm0[i], vih.ym0[i])
				//	if d < minD {
				//		minD = d
				//		minJ = j
				//	}
				//}
				//if minD < R {
				//	xc, yc := xlc[minJ], ylc[minJ]
				//	x0, y0 := vih.preXM0[i], vih.preYM0[i]
				//	x1, y1 := vih.xm0[i], vih.ym0[i]
				//
				//	a, b, c := x1-x0, y1-y0, (x0-x1)*xc+(y0-y1)*yc
				//	k := 1. / math.Sqrt(a*a+b*b)
				//
				//	d0 := k * (a*x0 + b*y0 + c)
				//	d1 := k * (a*x1 + b*y1 + c)
				//
				//	if d0*d1 < 0 {
				//		x0, y0 = x1, y1
				//	}
				//	dx := x0 - xc
				//	dy := y0 - yc
				//
				//	n := math.Sqrt(dx*dx + dy*dy)
				//	vih.xm0[i], vih.ym0[i] = xc+dx*R/n, yc+dy*R/n
				//}
			}
		}

		movingTime += time.Now().Sub(t0)
	}
	totalTime := time.Now().Sub(tt0)

	log.Printf("matrixTime:   %5.2f%%  %v\n", float64(matrixTime*100)/float64(totalTime), matrixTime)
	log.Printf("graphicsTime: %5.2f%%  %v\n", float64(graphicsTime*100)/float64(totalTime), graphicsTime)
	log.Printf("pointsTime:   %5.2f%%  %v\n", float64(pointsTime*100)/float64(totalTime), pointsTime)
	log.Printf("klyaksaTime:  %5.2f%%  %v\n", float64(klyaksaTime*100)/float64(totalTime), klyaksaTime)
	log.Printf("splinesTime:  %5.2f%%  %v\n", float64(splinesTime*100)/float64(totalTime), splinesTime)
	log.Printf("movingTime:   %5.2f%%  %v\n", float64(movingTime*100)/float64(totalTime), movingTime)
	log.Println("--------------------------")
	log.Printf("totalTime: %v\n", totalTime)

	f, _ := os.Create("img.gif")
	if err := gif.EncodeAll(f, &img); err != nil {
		log.Print(err)
	}
}

const (
	sizeX = 2400 * 30 / 44
	sizeY = 400
	k     = 200 / 2
	//minX, maxX = -2.5, 2.5
	//minY, maxY = -2.5, 2.5
	minX, maxX = -8., 22.
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

func (kl *Klyaksa) BuildSplinesX() (a, b, c, d []float64) {
	return kl.buildSplines(kl.x)
}

func (kl *Klyaksa) BuildSplinesY() (a, b, c, d []float64) {
	return kl.buildSplines(kl.y)
}

func (kl *Klyaksa) buildSplines(x []float64) (a, b, c, d []float64) {
	n := len(x)
	a = make([]float64, n)
	copy(a, x)

	c = kl.solve(x)

	//s := mat64.NewDense(n, n, nil)
	//s0 := mat64.NewDense(n, 1, nil)
	//for i := 1; i < n-1; i++ {
	//	s.Set(i, i-1, 1)
	//	s.Set(i, i, 4)
	//	s.Set(i, i+1, 1)
	//	s0.Set(i, 0, 3*((a[i+1]-a[i])-(a[i]-a[i-1])))
	//}
	//s.Set(0, n-1, 1)
	//s.Set(0, 0, 4)
	//s.Set(0, 1, 1)
	//s0.Set(0, 0, 3*((a[1]-a[0])-(a[0]-a[n-1])))
	//s.Set(n-1, n-2, 1)
	//s.Set(n-1, n-1, 4)
	//s.Set(n-1, 0, 1)
	//s0.Set(n-1, 0, 3*((a[0]-a[n-1])-(a[n-1]-a[n-2])))
	//
	//cc := &mat64.Dense{}
	//
	//if err := cc.Solve(s, s0); err != nil {
	//	panic(err)
	//}
	//c = cc.RawMatrix().Data[0:n]

	d = make([]float64, n)
	for i := 1; i < n; i++ {
		d[i] = (c[i] - c[i-1]) / 3
	}
	d[0] = (c[0] - c[n-1]) / 3

	b = make([]float64, n)
	for i := 1; i < n; i++ {
		b[i] = (a[i] - a[i-1]) + (2*c[i]+c[i-1])/3
	}
	b[0] = (a[0] - a[n-1]) + (2*c[0]+c[n-1])/3

	return
}

func (kl *Klyaksa) solveOld(x []float64) []float64 {
	n := len(x)

	s := mat64.NewDense(n, n, nil)
	s0 := mat64.NewDense(n, 1, nil)
	for i := 1; i < n-1; i++ {
		s.Set(i, i-1, 1)
		s.Set(i, i, 4)
		s.Set(i, i+1, 1)
		s0.Set(i, 0, 3*((x[i+1]-x[i])-(x[i]-x[i-1])))
	}
	s.Set(0, n-1, 1)
	s.Set(0, 0, 4)
	s.Set(0, 1, 1)
	s0.Set(0, 0, 3*((x[1]-x[0])-(x[0]-x[n-1])))
	s.Set(n-1, n-2, 1)
	s.Set(n-1, n-1, 4)
	s.Set(n-1, 0, 1)
	s0.Set(n-1, 0, 3*((x[0]-x[n-1])-(x[n-1]-x[n-2])))

	cc := &mat64.Dense{}

	if err := cc.Solve(s, s0); err != nil {
		panic(err)
	}
	return cc.RawMatrix().Data[0:n]
}

func (kl *Klyaksa) solve(x []float64) []float64 {
	n := len(x)

	a := make([]float64, n)
	b := make([]float64, n)
	c := make([]float64, n)
	f := make([]float64, n)

	for i := 1; i < n-1; i++ {
		a[i] = 1
		b[i] = 4
		c[i] = 1
		f[i] = 3 * ((x[i+1] - x[i]) - (x[i] - x[i-1]))
	}
	//a[0] = 0
	//c[n-1] = 0
	b[0], c[0], f[0] = 4, 1, 3*((x[1]-x[0])-(x[0]-x[n-1]))
	a[n-1], b[n-1], f[n-1] = 1, 4, 3*((x[0]-x[n-1])-(x[n-1]-x[n-2]))

	for i := 0; i < n-1; i++ {
		b[i+1] -= c[i] * a[i+1] / b[i]
		f[i+1] -= f[i] * a[i+1] / b[i]
		a[i+1] = 0
	}
	f[n-1] /= b[n-1]
	b[n-1] = 1
	for i := n - 2; i >= 0; i-- {
		f[i] -= f[i+1] * c[i]
		f[i] /= b[i]
		b[i] = 1
	}

	return f
}

func (kl *Klyaksa) BuildSplines() {
	kl.ax, kl.bx, kl.cx, kl.dx = kl.BuildSplinesX()
	kl.ay, kl.by, kl.cy, kl.dy = kl.BuildSplinesY()
}

func (kl *Klyaksa) AddSplineToImage(pal *image.Paletted) {
	ax, bx, cx, dx := kl.ax, kl.bx, kl.cx, kl.dx
	ay, by, cy, dy := kl.ay, kl.by, kl.cy, kl.dy

	for i := 0; i < len(kl.x); i++ {
		for t := 0.; t < 1; t += 0.01 {
			x := ax[i] - bx[i]*t + cx[i]*t*t - dx[i]*t*t*t
			y := ay[i] - by[i]*t + cy[i]*t*t - dy[i]*t*t*t
			xx, yy := convCoord(x, y)
			pal.SetColorIndex(xx, yy, 2)
		}
		xx, yy := convCoord(kl.x[i], kl.y[i])
		Tic(pal, xx, yy, 2)
		//pal.SetColorIndex(xx, yy, 1)
	}
}

func (kl *Klyaksa) AddToImage(pal *image.Paletted) {
	for i := 1; i < len(kl.x); i++ {
		x0, y0 := convCoord(kl.x[i-1], kl.y[i-1])
		x, y := convCoord(kl.x[i], kl.y[i])

		Line(pal, x0, y0, x, y, 2)
	}
	x0, y0 := convCoord(kl.x[0], kl.y[0])
	x, y := convCoord(kl.x[len(kl.x)-1], kl.y[len(kl.y)-1])

	Line(pal, x0, y0, x, y, 2)
}

func (kl *Klyaksa) AddPoints(maxDist float64) {
	ax, bx, cx, dx := kl.ax, kl.bx, kl.cx, kl.dx
	ay, by, cy, dy := kl.ay, kl.by, kl.cy, kl.dy

	var di int
	for i := 1; i < len(kl.x); i++ {
		x0, y0 := kl.x[i-1], kl.y[i-1]
		x, y := kl.x[i], kl.y[i]
		d := dist(x0, y0, x, y)
		n := int(d/maxDist) - 1
		if n <= 0 {
			continue
		}

		//fmt.Println(" -> ", x0, y0)
		var x1, y1, u1, v1 []float64
		for j := 0; j < n; j++ {
			t := float64(n-j) / float64(n+1)
			x := ax[i-di] - bx[i-di]*t + cx[i-di]*t*t - dx[i-di]*t*t*t
			y := ay[i-di] - by[i-di]*t + cy[i-di]*t*t - dy[i-di]*t*t*t
			//fmt.Println(t, x, y)
			x1 = append(x1, x)
			y1 = append(y1, y)

			//x1 = append(x1, (x0*float64(n-j)+x*float64(j+1))/float64(n+1))
			//y1 = append(y1, (y0*float64(n-j)+y*float64(j+1))/float64(n+1))
			u1 = append(u1, 0)
			v1 = append(v1, 0)
		}
		//fmt.Println(" <- ", x, y)

		x1 = append(x1, kl.x[i:]...)
		y1 = append(y1, kl.y[i:]...)
		u1 = append(u1, kl.u[i:]...)
		v1 = append(v1, kl.v[i:]...)

		kl.x = append(kl.x[:i], x1...)
		kl.y = append(kl.y[:i], y1...)
		kl.u = append(kl.u[:i], u1...)
		kl.v = append(kl.v[:i], v1...)
		i += n
		di += n
		//di += len(x1)
	}

	nn := len(kl.x)
	x0, y0 := kl.x[nn-1], kl.y[nn-1]
	x, y := kl.x[0], kl.y[0]
	d := dist(x0, y0, x, y)
	n := int(d/maxDist) - 1
	if n <= 0 {
		return
	}

	for j := 0; j < n; j++ {
		kl.x = append(kl.x, (x0*float64(n-j)+x*float64(j+1))/float64(n+1))
		kl.y = append(kl.y, (y0*float64(n-j)+y*float64(j+1))/float64(n+1))
		kl.u = append(kl.u, 0)
		kl.v = append(kl.v, 0)
	}
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

func addLabel(img *image.Paletted, x, y int, label string) {
	col := color.RGBA{200, 100, 0, 255}
	point := fixed.Point26_6{fixed.Int26_6(x * 64), fixed.Int26_6(y * 64)}

	d := &font.Drawer{
		Dst:  img,
		Src:  image.NewUniform(col),
		Face: basicfont.Face7x13,
		Dot:  point,
	}
	d.DrawString(label)
}
