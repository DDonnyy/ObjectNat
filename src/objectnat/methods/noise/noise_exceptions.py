class InvalidStepError(ValueError):
    def __init__(self, source_noise_db, target_noise_db, db_sim_step, div_, *args):
        if args:
            self.message = args[0]
        else:
            self.message = (
                f"The difference between `source_noise_db`({source_noise_db}) and `target_noise_db`({target_noise_db})"
                f" is not divisible by the step size ({db_sim_step}, remainder = {div_})"
            )

    def __str__(self):
        if self.message:
            return self.message
        return "The difference between `source_noise_db` and `target_noise_db` is not divisible by the step size"
