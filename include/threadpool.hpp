#ifndef threadpool_h
#define threadpool_h

#include <atomic>
#include <condition_variable>
#include <functional>
#include <future>
#include <mutex>
#include <queue>
#include <thread>
#include <typeinfo>
#include <type_traits>
#include <vector>

/* ThreadPool class */
class ThreadPool
{
public:
    // Constructor
    ThreadPool(int8_t numThreads)
    {
        createThreads(numThreads);
    }

    // No copy c'tor and assignment
    ThreadPool& operator = (const ThreadPool&) = delete;
    ThreadPool(const ThreadPool&) = delete;

    // Destructor
    ~ThreadPool()
    {
        shutdown = true;
        notifier.notify_all();
        for (size_t i = 0; i < threads.size(); ++i)
        {
            threads[i].join();
        }
    }

    //add any arg # function to queue
    template <typename Func, typename... Args>
    auto push(Func &&f, Args &&... args)
    {
        //get return type of the function
        using RetType = std::invoke_result_t<decltype(f),Args...>;
        std::packaged_task<RetType()> task(std::bind(f, std::forward<Args>(args)...));

        std::future<RetType> future = task.get_future();

        {
            // lock jobQueue mutex, add job to the job queue
            std::unique_lock<std::mutex> lock(JobMutex);

            //place the job into the queue
            jobQueue.emplace(std::packaged_task<void()>(std::move(task)));
        }
        notifier.notify_one();
        return future;
    }

private:

    std::vector<std::thread> threads;
    std::queue<std::packaged_task<void()>> jobQueue;
    std::condition_variable notifier;
    std::mutex JobMutex;
    std::atomic<bool> shutdown = false;

    void createThreads(uint8_t numThreads)
    {
        auto threadFunc = [this]() {
            while (true)
            {
                std::packaged_task<void()> job;

                {
                    std::unique_lock<std::mutex> lock(JobMutex);
                    notifier.wait(lock, [this] { return !jobQueue.empty() || shutdown; });

                    if (shutdown)
                    {
                        break;
                    }

                    job = std::move(jobQueue.front());

                    jobQueue.pop();
                }
                job();
            }
        };

        threads.reserve(numThreads);
        for (int i = 0; i != numThreads; ++i)
        {
            threads.emplace_back(threadFunc);
        }
    }

}; /* end ThreadPool Class */

#endif // threadpool_h