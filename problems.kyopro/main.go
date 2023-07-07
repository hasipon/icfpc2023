package main

import (
	"encoding/json"
	"fmt"
	"log"
	"os"
)

type Attendee struct {
	X float64 `json:"x"`
	Y float64 `json:"y"`
}

type Problem struct {
	RoomWidth       float64     `json:"room_width"`
	RoomHeight      float64     `json:"room_height"`
	StageWidth      float64     `json:"stage_width"`
	StageHeight     float64     `json:"stage_height"`
	StageBottomLeft []float64   `json:"stage_bottom_left"`
	Musicians       []int       `json:"musicians"`
	Attendees       []*Attendee `json:"attendees"`
}

// func (p *Problem) Stringer() string {}
func (p *Problem) show(file *os.File) {
	fmt.Fprintf(file, "%f %f\n", p.RoomWidth, p.RoomHeight)
	fmt.Fprintf(file, "%f %f\n", p.StageWidth, p.StageHeight)
	fmt.Fprintf(file, "%f %f\n", p.StageBottomLeft[0], p.StageBottomLeft[1])
	fmt.Fprintf(file, "%d\n", len(p.Musicians))
	for _, m := range p.Musicians {
		fmt.Fprintf(file, "%d ", m)
	}
	fmt.Fprintf(file, "\n")
	fmt.Fprintf(file, "%d\n", len(p.Attendees))
	for _, a := range p.Attendees {
		fmt.Fprintf(file, "%f %f\n", a.X, a.Y)
	}
}

func main() {
	for i := 1; i <= 45; i++ {
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
