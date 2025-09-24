from openai import OpenAI

# Set OpenAI's API key and API base to use vLLM's API server.
openai_api_key = "EMPTY"
openai_api_base = "http://localhost:8000/v1"

# Define the prompts to send
prompts = [
    {"role": "system", "content": "You are a helpful assistant."},
    {"role": "user", "content": "Tell me a joke."},
]

# Initialize the client
client = OpenAI(
    api_key=openai_api_key,
    base_url=openai_api_base,
)

# Send the prompts and print the responses
chat_response = client.chat.completions.create(
    model="meta-llama/Llama-2-7b-chat-hf",
    messages=prompts
)
reply = chat_response.choices[0].message.content
for prompt in prompts:
    role = prompt["role"].capitalize()
    print(f"{role}: {prompt['content']}")
print(f"Assistant: {reply}")

