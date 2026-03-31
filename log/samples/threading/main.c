#include <pthread.h>
#include <signal.h>
#include <unistd.h>
#include <stdio.h>

#include <log.h>

static bool s_stop = false;
pthread_mutex_t MUTEX_LOG;

// Mutex lock function for log module
void log_lock(bool lock, void* mutex) {
    pthread_mutex_t *LOCK = (pthread_mutex_t*)(mutex);

    if (lock)
    {
	pthread_mutex_lock(LOCK);
    }
    else
    {
	pthread_mutex_unlock(LOCK);
    }
}

// Thread runner function
void* run(void *arg) {
    // Retrieve this thread's id
    int id = *(int*)arg;

    log_debug("+++++id: %d", id);
    unsigned int idx = 0;

    // Forever
    for (;;) {
	// Log the thread's id and the iteration index
	log_debug("thread id: %d idx:%u", id, idx++);

	// Sleep 100ms
	    usleep(1000 * 100);

	    if (s_stop) { break; }
    }

    return 0;
}

// Signal catcher
void stop(int signo) {
    (void)signo;
    s_stop = true;
}

int main() {
    // Create a mutex
    pthread_mutex_init(&MUTEX_LOG, NULL);

    // Set the lock function and the mutex for the log library
    log_set_lock(log_lock, &MUTEX_LOG);

    FILE* fp = fopen ("./demo.log", "a+");

    log_set_level(LOGC_DEBUG);
    log_add_fp(fp, LOGC_DEBUG);

    // Redirect signals to stop() fucntion
    signal(SIGINT, stop);
    signal(SIGKILL, stop);
    signal(SIGTERM, stop);

    pthread_t thread_handles[8];
    int id_array[8] = {0,1,2,3,4,5,6,7};

    // Create 8 threads, each passing their id as parameter
    for (int i= 0; i<8; i++) {
	pthread_create(&thread_handles[i], NULL, run, &id_array[i]);
    }

    void* thread_retvals[8] = {NULL};

    // Start threads
    for (int i=0; i<8; i++) {
	pthread_join(thread_handles[i], &thread_retvals[i]);
    }

    pthread_mutex_destroy(&MUTEX_LOG);

    fclose(fp);

    return 0;
}
