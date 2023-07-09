package main

import (
	"context"
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"log"
	"os"
	"os/exec"
	"path"
	"path/filepath"
	"strconv"
	"strings"
	"sync"
	"sync/atomic"
	"time"
)

var (
	RepoRoot          string
	OutputDir         string
	ActionID          string
	SolverName        string
	SolverPath        string
	TimeoutEntire     time.Duration
	TimeoutPerProblem time.Duration
	MaxParallel       int
	ProblemIDs        []int
	AutoCommit        bool

	runsMtx    sync.RWMutex
	solverRuns []*SolverRun
)

type SolverRun struct {
	ProblemID        int
	SolverPath       string
	Args             []string
	Timeout          time.Duration
	StartAt          time.Time
	SolutionFileName string
	OnFinish         func(*SolverRun)

	Running    atomic.Bool
	StdOutPath string
	StdErrPath string
	Error      error
	ExitCode   int
}

func (run *SolverRun) Run(ctx context.Context) {
	if run.OnFinish != nil {
		defer run.OnFinish(run)
	}

	ctx, cancel := context.WithTimeout(ctx, run.Timeout)
	defer cancel()

	run.Running.Store(true)
	defer func() {
		run.Running.Store(false)
	}()

	cmd := exec.CommandContext(ctx, run.SolverPath, run.Args...)
	cmd.Env = append(
		os.Environ(),
		fmt.Sprintf("PROBLEM_ID=%v", run.ProblemID),
		fmt.Sprintf("REPO_ROOT=%v", RepoRoot),
		fmt.Sprintf("TIMEOUT=%d", int(run.Timeout.Seconds())),
	)

	cmdStr := run.SolverPath + " " + strings.Join(run.Args, " ")
	log.Println("Run:", cmdStr, "ProblemID:", run.ProblemID)

	run.Error = func() error {
		problem, err := os.ReadFile(path.Join(RepoRoot, "problems.kyopro", fmt.Sprintf("%d.kyopro", run.ProblemID)))
		if err != nil {
			return errors.New(fmt.Sprintf("execCommand os.ReadFile Error: %v", err))
		}

		stdoutReader, err := cmd.StdoutPipe()
		if err != nil {
			return errors.New(fmt.Sprintf("execCommand StdoutPipe Error: %v", err))
		}
		stderrReader, err := cmd.StderrPipe()
		if err != nil {
			return errors.New(fmt.Sprintf("execCommand StderrPipe Error: %v", err))
		}
		stdinWriter, err := cmd.StdinPipe()
		if err != nil {
			return errors.New(fmt.Sprintf("execCommand StdinPipe Error: %v", err))
		}
		defer stdinWriter.Close()

		stdoutFilePath := path.Join(OutputDir, fmt.Sprintf("%d-%s-stdout.txt", run.ProblemID, filepath.Base(run.SolverPath)))
		stdoutFile, err := os.OpenFile(stdoutFilePath, os.O_RDWR|os.O_CREATE, 0644)
		if err != nil {
			return errors.New(fmt.Sprintf("execCommand OpenFile Error: %v", err))
		}
		defer stdoutFile.Close()
		run.StdOutPath = stdoutFilePath

		stderrFilePath := path.Join(OutputDir, fmt.Sprintf("%d-%s-stderr.txt", run.ProblemID, filepath.Base(run.SolverPath)))
		stderrFile, err := os.OpenFile(stderrFilePath, os.O_RDWR|os.O_CREATE, 0644)
		if err != nil {
			return errors.New(fmt.Sprintf("execCommand OpenFile Error: %v", err))
		}
		defer stderrFile.Close()
		run.StdErrPath = stderrFilePath

		ch := make(chan error)
		go func() {
			ch <- func() error {
				_, err := io.Copy(stdoutFile, stdoutReader)
				if err != nil {
					return errors.New(fmt.Sprintf("[stdout io.Copy error] %v", err))
				}
				return nil
			}()
		}()

		go func() {
			ch <- func() error {
				_, err := io.Copy(stderrFile, stderrReader)
				if err != nil {
					return errors.New(fmt.Sprintf("[stdout io.Copy error] %v", err))
				}
				return nil
			}()
		}()

		err = cmd.Start()
		if err != nil {
			return errors.New(fmt.Sprintf("execCommand Start Error: %v", err))
		}

		_, err = stdinWriter.Write(problem)
		if err != nil {
			return errors.New(fmt.Sprintf("execCommand stdinWrite Error: %v", err))
		}

		err1 := <-ch
		err2 := <-ch
		if err1 != nil {
			if err2 != nil {
				return errors.New(fmt.Sprintf("%v %v", err1, err2))
			} else {
				return err1
			}
		} else if err2 != nil {
			return err2
		}

		err = cmd.Wait()
		run.ExitCode = cmd.ProcessState.ExitCode()
		if err != nil {
			return errors.New(fmt.Sprintf("execCommand Wait Error: %v", err))
		}
		return nil
	}()
}

func gitAdd(newFileName string) {
	cmd := exec.Command("git", "add", "solutions/"+newFileName)
	cmd.Dir = RepoRoot
	if out, err := cmd.CombinedOutput(); err != nil {
		log.Println("[gitAdd] exec.Command Error", err)
		return
	} else {
		log.Println(string(out))
	}
}

func gitCommitAndPush() {
	commitMessage := fmt.Sprintf("Runner Action:%v Add solution of %v", ActionID, SolverName)
	cmd := exec.Command("git", "commit", "-m", fmt.Sprintf("%q", commitMessage))
	cmd.Dir = RepoRoot
	if out, err := cmd.CombinedOutput(); err != nil {
		log.Println("[gitCommitAndPush] exec.Command Error", err)
		return
	} else {
		log.Println(string(out))
	}

	cmd = exec.Command("git", "push")
	cmd.Dir = RepoRoot
	if out, err := cmd.CombinedOutput(); err != nil {
		log.Println("[gitCommitAndPush] exec.Command Error", err)
		return
	} else {
		log.Println(string(out))
	}
}

func handleWorkerEnd(run *SolverRun) {
	runsMtx.RLock()
	defer runsMtx.RUnlock()

	// copy solution
	if run.StdOutPath != "" {
		stdout, err := os.ReadFile(run.StdOutPath)
		if err != nil {
			log.Println("os.ReadFile(run.StdOutPath) Error:", err)
		} else if !json.Valid(stdout) {
			log.Println("Invalid JSON output:", string(stdout))
		} else {
			filePath := path.Join(RepoRoot, "solutions", run.SolutionFileName)
			err := os.WriteFile(filePath, stdout, 0644)
			if err != nil {
				log.Println("os.WriteFile Error:", err, filePath)
			} else {
				log.Println("Saved:", filePath)
				if AutoCommit {
					gitAdd(run.SolutionFileName)
				}
			}
		}
	}
}

func fileExists(filename string) bool {
	info, err := os.Stat(filename)
	if os.IsNotExist(err) {
		return false
	}
	return !info.IsDir()
}

func mainRun() {
	ctx, cancel := context.WithTimeout(context.Background(), TimeoutEntire)
	defer cancel()

	runsMtx.Lock()
	for _, problemID := range ProblemIDs {
		actionSuffix := ""
		if ActionID != "" {
			actionSuffix = "_" + ActionID
		}
		solutionFileName := fmt.Sprintf("%d-%s%s.json", problemID, SolverName, actionSuffix)
		solutionFilePath := path.Join(RepoRoot, "solutions", solutionFileName)
		if fileExists(solutionFilePath) {
			log.Println(solutionFileName, "already exists. skip execution.")
			continue
		}

		run := &SolverRun{
			ProblemID:        problemID,
			SolverPath:       SolverPath,
			Args:             nil,
			Timeout:          TimeoutPerProblem,
			SolutionFileName: solutionFileName,
			OnFinish:         handleWorkerEnd,
		}
		solverRuns = append(solverRuns, run)
	}
	runsMtx.Unlock()

	sema := make(chan bool, MaxParallel)
	for i := 0; i < MaxParallel; i++ {
		sema <- true
	}

	ticker := time.NewTicker(30 * time.Second)
	wg := &sync.WaitGroup{}
	nextIndex := 0
	for {
		select {
		case <-ticker.C:
			if AutoCommit {
				gitCommitAndPush()
			}

		case <-sema:
			runsMtx.RLock()
			if nextIndex < len(solverRuns) {
				wg.Add(1)
				go func(r *SolverRun) {
					defer func() {
						sema <- true
						wg.Done()
					}()
					r.Run(ctx)
				}(solverRuns[nextIndex])
				nextIndex++
				log.Printf("Enqueue run [%d/%d]\n", nextIndex, len(solverRuns))
			} else {
				cancel()
			}
			runsMtx.RUnlock()

		case <-ctx.Done():
			if ctx.Err() != nil {
				log.Println(ctx.Err())
			}
			log.Println("Waiting all solver to stop")
			wg.Wait()
			if AutoCommit {
				gitCommitAndPush()
			}
			log.Println("Complete")
			return
		}
	}

}

func main() {
	log.Println("Runner Started")

	RepoRoot = os.Getenv("REPO_ROOT")
	if RepoRoot == "" {
		fmt.Println("env REPO_ROOT must be set")
		os.Exit(1)
	}

	SolverName = os.Getenv("SOLVER_NAME")
	if SolverName == "" {
		fmt.Println("env SOLVER_NAME must be set")
		os.Exit(1)
	}

	SolverPath = os.Getenv("SOLVER_PATH")
	if SolverPath == "" {
		fmt.Println("env SOLVER_PATH must be set")
		os.Exit(1)
	}

	ActionID = os.Getenv("ACTION_ID")
	OutputDir = os.Getenv("OUTPUT_DIR")
	if OutputDir == "" {
		OutputDir = "./output"
	}

	TimeoutEntire = time.Hour
	if os.Getenv("TIMEOUT_ENTIRE") != "" {
		d, err := time.ParseDuration(os.Getenv("TIMEOUT_ENTIRE"))
		if err != nil {
			log.Fatal("TIMEOUT_ENTIRE", err)
		}
		TimeoutEntire = d
	}

	TimeoutPerProblem = time.Minute
	if os.Getenv("TIMEOUT_PER_PROBLEM") != "" {
		d, err := time.ParseDuration(os.Getenv("TIMEOUT_PER_PROBLEM"))
		if err != nil {
			log.Fatal("TIMEOUT_PER_PROBLEM", err)
		}
		TimeoutPerProblem = d
	}

	MaxParallel = 4
	if os.Getenv("MAX_P") != "" {
		n, err := strconv.Atoi(os.Getenv("MAX_P"))
		if err != nil {
			log.Fatal("MAX_P", err)
		}
		if 0 < n && n <= 20 {
			MaxParallel = n
		}
	}

	targetIDs := "1-10"
	if os.Getenv("PROBLEM_IDS") != "" {
		targetIDs = os.Getenv("PROBLEM_IDS")
	}

	s := strings.TrimSpace(targetIDs)
	for _, e := range strings.Split(s, " ") {
		if strings.Contains(e, "-") {
			ab := strings.SplitN(s, "-", 2)
			a, err := strconv.Atoi(ab[0])
			if err != nil {
				log.Fatal("PROBLEM_IDS a-b", err, s)
			}
			b, err := strconv.Atoi(ab[1])
			if err != nil {
				log.Fatal("PROBLEM_IDS a-b", err, s)
			}
			if a > b {
				if err != nil {
					log.Fatal("PROBLEM_IDS a>b", err, s)
				}
			}

			for i := a; i <= b; i++ {
				ProblemIDs = append(ProblemIDs, i)
			}
		} else {
			a, err := strconv.Atoi(e)
			if err != nil {
				log.Fatal("PROBLEM_IDS a", err, s)
			}
			ProblemIDs = append(ProblemIDs, a)
		}
	}

	AutoCommit = false
	if os.Getenv("AUTO_COMMIT") != "" {
		AutoCommit = os.Getenv("AUTO_COMMIT") != "0" && os.Getenv("AUTO_COMMIT") != "no"
	}

	log.Println("==========")
	log.Println("RepoRoot", RepoRoot)
	log.Println("ActionID", ActionID)
	log.Println("OutputDir", OutputDir)
	log.Println("MaxParallel", MaxParallel)
	log.Println("TimeoutEntire", TimeoutEntire)
	log.Println("TimeoutPerProblem", TimeoutPerProblem)
	log.Println("ProblemIDs", ProblemIDs)
	log.Println("AutoCommit", AutoCommit)
	log.Println("==========")

	err := os.MkdirAll(OutputDir, 0755)
	if err != nil {
		log.Fatalf("os.MkdirAll: %s", err)
	}

	mainRun()
}
