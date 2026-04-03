from openai import OpenAI
import os
import time

# Get the model from the environment variable
model = os.getenv("MODEL", "meta-llama/Llama-3.1-8B-Instruct")

# Set OpenAI's API key and API base to use vLLM's API server.
openai_api_key = "EMPTY"
openai_api_base = "http://localhost:8000/v1"

# Define the prompts to send
prompts = [
    {"role": "system", "content": "You are a helpful assistant."},
    {"role": "user", "content": "Introduce yourself."},
]

# Initialize the client
client = OpenAI(
    api_key=openai_api_key,
    base_url=openai_api_base,
)

# Wait for vLLM server to be ready
timeout = 300
interval = 10
start = time.time()

while time.time() - start < timeout:
    try:
        chat_response = client.chat.completions.create(
            model=model,
            messages=prompts
        )
        break
    except Exception as e:
        elapsed = int(time.time() - start)
        remaining = timeout - elapsed
        print(f"vLLM server not ready ({e}). Retrying in {interval}s... ({remaining}s remaining)")
        time.sleep(interval)
else:
    raise TimeoutError(f"vLLM server did not respond within {timeout}s. You can try again.")

reply = chat_response.choices[0].message.content
for prompt in prompts:
    role = prompt["role"].capitalize()
    print(f"{role}: {prompt['content']}")
print(f"Assistant: {reply}")

