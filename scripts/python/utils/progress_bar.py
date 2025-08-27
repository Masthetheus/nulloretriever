"""Manual implemented progress bar for global purposes"""
import datetime as dt

def progress_bar(current, total, start, size=30):
    """Progress bar for general completion tracking
    Args:
        current(int): time related to function call
        total(int): total of units to be processed (organisms, k values, etc)
        start(int): time at the first function call
        size(int): total bar size
    Returns:
        terminal: bar displaying estimated time of completion
    """
    now = time.time()
    elapsed = now - start
    if current > 0:
        estimated_time = elapsed / current
        remainder = estimated_time * (total - current)
    else:
        remainder = 0
    minutes, seconts = divmod(int(remainder), 60)
    proportion = current / total
    complete = int(proportion * size)
    bar = "|" + "o" * complete + "-" * (size - complete) + "|"
    print(f"\r{bar} {current}/{total} - Estimated: {minutes:02d}:{seconts:02d} left", end='', flush=True)
    if current == total:
        print()
