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

	"github.com/caarlos0/env/v9"
)

type Config struct {
	RepoRoot          string        `env:"REPO_ROOT"`
	OutputDir         string        `env:"OUTPUT_DIR" envDefault:"./output"`
	ActionID          string        `env:"ACTION_ID"`
	SolverName        string        `env:"SOLVER_NAME"`
	SolverPath        string        `env:"SOLVER_PATH"`
	ProblemIDRange    string        `env:"PROBLEM_IDS" envDefault:"1"`
	TimeoutEntire     time.Duration `env:"TIMEOUT_ENTIRE" envDefault:"1h"`
	TimeoutPerProblem time.Duration `env:"TIMEOUT_PER_PROBLEM" envDefault:"1m"`
	MaxParallel       int           `env:"MAX_PARALLEL" envDefault:"4"`
	AutoCommit        bool          `env:"AUTO_COMMIT" envDefault:"0"`
}

var (
	conf       Config
	ProblemIDs []int
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
		fmt.Sprintf("REPO_ROOT=%v", conf.RepoRoot),
		fmt.Sprintf("TIMEOUT=%d", int(run.Timeout.Seconds())),
	)

	cmdStr := run.SolverPath + " " + strings.Join(run.Args, " ")
	log.Println("Run:", cmdStr, "ProblemID:", run.ProblemID)

	run.Error = func() error {
		problem, err := os.ReadFile(path.Join(conf.RepoRoot, "problems.kyopro", fmt.Sprintf("%d.kyopro", run.ProblemID)))
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

		stdoutFilePath := path.Join(conf.OutputDir, fmt.Sprintf("%d-%s-stdout.txt", run.ProblemID, filepath.Base(run.SolverPath)))
		stdoutFile, err := os.OpenFile(stdoutFilePath, os.O_RDWR|os.O_CREATE, 0644)
		if err != nil {
			return errors.New(fmt.Sprintf("execCommand OpenFile Error: %v", err))
		}
		defer stdoutFile.Close()
		run.StdOutPath = stdoutFilePath

		stderrFilePath := path.Join(conf.OutputDir, fmt.Sprintf("%d-%s-stderr.txt", run.ProblemID, filepath.Base(run.SolverPath)))
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
	cmd.Dir = conf.RepoRoot
	if out, err := cmd.CombinedOutput(); err != nil {
		log.Println("[gitAdd] exec.Command Error", err)
		return
	} else {
		log.Println(string(out))
	}
}

func gitCommitAndPush() {
	commitMessage := fmt.Sprintf("Runner Action:%v Add solution of %v", conf.ActionID, conf.SolverName)
	cmd := exec.Command("git", "commit", "-m", fmt.Sprintf("%q", commitMessage))
	cmd.Dir = conf.RepoRoot
	if out, err := cmd.CombinedOutput(); err != nil {
		log.Println("[gitCommitAndPush] exec.Command Error", err)
		return
	} else {
		log.Println(string(out))
	}

	cmd = exec.Command("git", "push")
	cmd.Dir = conf.RepoRoot
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
			filePath := path.Join(conf.RepoRoot, "solutions", run.SolutionFileName)
			err := os.WriteFile(filePath, stdout, 0644)
			if err != nil {
				log.Println("os.WriteFile Error:", err, filePath)
			} else {
				log.Println("Saved:", filePath)
				if conf.AutoCommit {
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
	ctx, cancel := context.WithTimeout(context.Background(), conf.TimeoutEntire)
	defer cancel()

	runsMtx.Lock()
	for _, problemID := range ProblemIDs {
		actionSuffix := ""
		if conf.ActionID != "" {
			actionSuffix = "_" + conf.ActionID
		}
		solutionFileName := fmt.Sprintf("%d-%s%s.json", problemID, conf.SolverName, actionSuffix)
		solutionFilePath := path.Join(conf.RepoRoot, "solutions", solutionFileName)
		if fileExists(solutionFilePath) {
			log.Println(solutionFileName, "already exists. skip execution.")
			continue
		}

		run := &SolverRun{
			ProblemID:        problemID,
			SolverPath:       conf.SolverPath,
			Args:             nil,
			Timeout:          conf.TimeoutPerProblem,
			SolutionFileName: solutionFileName,
			OnFinish:         handleWorkerEnd,
		}
		solverRuns = append(solverRuns, run)
	}
	runsMtx.Unlock()

	sema := make(chan bool, conf.MaxParallel)
	for i := 0; i < conf.MaxParallel; i++ {
		sema <- true
	}

	ticker := time.NewTicker(30 * time.Second)
	wg := &sync.WaitGroup{}
	nextIndex := 0
	for {
		select {
		case <-ticker.C:
			if conf.AutoCommit {
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
			}
			last := nextIndex == len(solverRuns)
			runsMtx.RUnlock()

			if last {
				wg.Wait()
				cancel()
			}

		case <-ctx.Done():
			if ctx.Err() != nil {
				log.Println(ctx.Err())
			}
			log.Println("Waiting all solver to stop")
			wg.Wait()
			if conf.AutoCommit {
				gitCommitAndPush()
			}
			log.Println("Complete")
			return
		}
	}
}

func parseProblemIDRange() {
	for _, e := range strings.Split(conf.ProblemIDRange, " ") {
		if strings.Contains(e, "-") {
			ab := strings.SplitN(e, "-", 2)
			a, err := strconv.Atoi(ab[0])
			if err != nil {
				log.Fatal("PROBLEM_IDS a-b", err, e)
			}
			b, err := strconv.Atoi(ab[1])
			if err != nil {
				log.Fatal("PROBLEM_IDS a-b", err, e)
			}
			if a > b {
				if err != nil {
					log.Fatal("PROBLEM_IDS a>b", err, e)
				}
			}

			for i := a; i <= b; i++ {
				ProblemIDs = append(ProblemIDs, i)
			}
		} else {
			a, err := strconv.Atoi(e)
			if err != nil {
				log.Fatal("PROBLEM_IDS a", err, e)
			}
			ProblemIDs = append(ProblemIDs, a)
		}
	}
}

func main() {
	log.Println("Runner Started")

	if err := env.Parse(&conf); err != nil {
		log.Fatal("%+v\n", err)
	}
	parseProblemIDRange()

	log.Println("==========")
	log.Println("RepoRoot", conf.RepoRoot)
	log.Println("SolverName", conf.SolverName)
	log.Println("SolverPath", conf.SolverPath)
	log.Println("ActionID", conf.ActionID)
	log.Println("OutputDir", conf.OutputDir)
	log.Println("MaxParallel", conf.MaxParallel)
	log.Println("TimeoutEntire", conf.TimeoutEntire)
	log.Println("TimeoutPerProblem", conf.TimeoutPerProblem)
	log.Println("ProblemIDs", ProblemIDs)
	log.Println("AutoCommit", conf.AutoCommit)
	log.Println("==========")

	if conf.RepoRoot == "" {
		log.Fatal("RepoRoot is not set")
	}
	if conf.SolverPath == "" {
		log.Fatal("SolverPath is not set")
	}
	if !fileExists(conf.SolverPath) {
		if fileExists(path.Join(conf.RepoRoot, conf.SolverPath)) {
			conf.SolverPath = path.Join(conf.RepoRoot, conf.SolverPath)
		} else {
			log.Fatal("File not exist at SolverPath", conf.SolverPath)
		}
	}
	if conf.SolverName == "" {
		log.Fatal("SolverName is not set")
	}

	err := os.MkdirAll(conf.OutputDir, 0755)
	if err != nil {
		log.Fatalf("os.MkdirAll: %s", err)
	}

	mainRun()
}
