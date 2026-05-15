#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <string_view>

#include <gtest/gtest.h>

namespace
{
std::filesystem::path unique_temp_path(std::string_view stem,
				       std::string_view extension)
{
	return std::filesystem::temp_directory_path() /
	       (std::string(stem) + "_" +
		std::to_string(std::chrono::steady_clock::now()
				       .time_since_epoch()
				       .count()) +
		std::string(extension));
}

std::filesystem::path write_gfa_fixture(std::string_view stem,
					std::string_view contents)
{
	std::filesystem::path fp = unique_temp_path(stem, ".gfa");
	std::ofstream out(fp);
	out << contents;
	return fp;
}

std::string read_file(const std::filesystem::path &fp)
{
	std::ifstream in(fp);
	return std::string(std::istreambuf_iterator<char>(in),
			   std::istreambuf_iterator<char>());
}

struct PovuRun {
	int status;
	std::string stderr_text;
};

PovuRun run_gfa2vcf(const std::filesystem::path &gfa_fp)
{
	const std::filesystem::path povu =
		std::filesystem::path(POVU_SOURCE_DIR) / "bin" / "povu";
	const std::filesystem::path stderr_fp =
		unique_temp_path("povu_gfa2vcf_stderr", ".txt");

	pid_t pid = fork();
	if (pid == 0) {
		int stdout_fd = open("/dev/null", O_WRONLY);
		int stderr_fd = open(stderr_fp.c_str(),
				     O_WRONLY | O_CREAT | O_TRUNC, 0600);
		if (stdout_fd < 0 || stderr_fd < 0)
			_exit(126);

		dup2(stdout_fd, STDOUT_FILENO);
		dup2(stderr_fd, STDERR_FILENO);
		close(stdout_fd);
		close(stderr_fd);

		execl(povu.c_str(), povu.c_str(), "gfa2vcf", "-i",
		      gfa_fp.c_str(), "-t", "1", "-P", "HG1",
		      static_cast<char *>(nullptr));
		_exit(127);
	}

	if (pid < 0) {
		return PovuRun{127, "fork failed"};
	}

	int status = 0;
	waitpid(pid, &status, 0);

	std::string stderr_text = read_file(stderr_fp);
	std::filesystem::remove(stderr_fp);
	return PovuRun{status, stderr_text};
}

void expect_controlled_rejection(const std::filesystem::path &gfa_fp)
{
	PovuRun run = run_gfa2vcf(gfa_fp);
	EXPECT_TRUE(WIFEXITED(run.status)) << run.stderr_text;
	if (WIFEXITED(run.status))
		EXPECT_NE(WEXITSTATUS(run.status), 0) << run.stderr_text;
	EXPECT_FALSE(WIFSIGNALED(run.status))
		<< "signal " << WTERMSIG(run.status) << "\n"
		<< run.stderr_text;
	EXPECT_NE(run.stderr_text.find("Invalid GFA"), std::string::npos)
		<< run.stderr_text;
}
} // namespace

TEST(Gfa2VcfTest, RejectsUnsupportedRecordWithoutSignal)
{
	const std::filesystem::path gfa_fp = write_gfa_fixture(
		"povu_unsupported_record",
		"H\tVN:Z:1.0\n"
		"Z\tnot\tvalid\n"
		"S\t0\tA\n"
		"P\tHG1#1#chr1\t0+\t*\n");

	expect_controlled_rejection(gfa_fp);
	std::filesystem::remove(gfa_fp);
}

TEST(Gfa2VcfTest, RejectsSegmentMissingSequenceWithoutSignal)
{
	const std::filesystem::path gfa_fp = write_gfa_fixture(
		"povu_missing_sequence",
		"H\tVN:Z:1.0\n"
		"S\t0\n"
		"P\tHG1#1#chr1\t0+\t*\n");

	expect_controlled_rejection(gfa_fp);
	std::filesystem::remove(gfa_fp);
}
