package main

import (
	"encoding/json"
	"fmt"
	"log"
	"os"
)

type Attendee struct {
	X      float64   `json:"x"`
	Y      float64   `json:"y"`
	Teasts []float64 `json:"tastes"`
}

type Pillar struct {
	Center []float64 `json:"center"`
	Radius float64   `json:"radius"`
}

type Problem struct {
	RoomWidth       float64     `json:"room_width"`
	RoomHeight      float64     `json:"room_height"`
	StageWidth      float64     `json:"stage_width"`
	StageHeight     float64     `json:"stage_height"`
	StageBottomLeft []float64   `json:"stage_bottom_left"`
	Musicians       []int       `json:"musicians"`
	Attendees       []*Attendee `json:"attendees"`
	Pillars         []*Pillar   `json:"pillars"`
}

// func (p *Problem) Stringer() string {}
func (p *Problem) show(file *os.File) {
	fmt.Fprintf(file, "%f %f\n", p.RoomWidth, p.RoomHeight)
	fmt.Fprintf(file, "%f %f\n", p.StageWidth, p.StageHeight)
	fmt.Fprintf(file, "%f %f\n", p.StageBottomLeft[0], p.StageBottomLeft[1])
	cnt := map[int]bool{}
	for _, m := range p.Musicians {
		cnt[m] = true
	}
	fmt.Fprintf(file, "%d %d\n", len(p.Musicians), len(cnt))
	for _, m := range p.Musicians {
		fmt.Fprintf(file, "%d ", m)
		cnt[m] = true
	}
	fmt.Fprintf(file, "\n")
	fmt.Fprintf(file, "%d\n", len(p.Attendees))
	for _, a := range p.Attendees {
		fmt.Fprintf(file, "%f %f", a.X, a.Y)
		for _, t := range a.Teasts {
			fmt.Fprintf(file, " %f", t)
		}
		fmt.Fprintf(file, "\n")
	}
	fmt.Fprintf(file, "%d\n", len(p.Pillars))
	for _, a := range p.Pillars {
		fmt.Fprintf(file, "%f %f %f\n", a.Center[0], a.Center[1], a.Radius)
	}
}

func main() {
	for i := 1; i <= 90; i++ {
		fmt.Println(i)
		initJson, err := os.ReadFile(fmt.Sprintf("../problem.json/%d.json", i))
		if err != nil {
			panic(err)
			log.Fatalln(err)
		}
		outFile, err := os.Create(fmt.Sprintf("../problems.kyopro/%d.kyopro", i))
		problem := Problem{}
		err = json.Unmarshal(initJson, &problem)
		if err != nil {
			log.Fatalln(err)
		}
		problem.show(outFile)
	}
}
